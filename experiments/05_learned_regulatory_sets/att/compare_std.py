import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import csv
import glob
import json
import os
import sys
import numpy as np
from scipy.stats import norm

sys.path.insert(0, _config_path('${PROJECT_ROOT}/work/fset/att'))
import att_common_std as C

NFOLD = 5
COLS = [('RAW_ALL', C.C_RAW_ALL), ('STD_ALL', C.C_STD_ALL),
        ('STD_COMMON', C.C_STD_COMMON), ('STD_LOW', C.C_STD_LOW),
        ('STD_RARE', C.C_STD_RARE), ('STD_VRARE', C.C_STD_VRARE),
        ('RAW_COMMON', C.C_RAW_COMMON), ('STD_CCRE', C.C_STD_CCRE),
        ('STD_NOCCRE', C.C_STD_NOCCRE)]

def load_pheno():
    from cyvcf2 import VCF
    v = VCF(C.VCFD + '/chr22.vcf.gz')
    pos = dict((s, i) for i, s in enumerate(v.samples))
    idx, ys, cs = [], [], []
    for r in csv.DictReader(open(C.PHENO), delimiter='\t'):
        i = pos.get(r['sample_id'])
        if i is None:
            continue
        try:
            yv = float(r['y']); cv = [float(r[c]) for c in C.COVARS]
        except (ValueError, TypeError, KeyError):
            continue
        if not np.isfinite(yv) or not np.isfinite(cv).all():
            continue
        idx.append(i); ys.append(yv); cs.append(cv)
    idx = np.asarray(idx, dtype=np.int64); o = np.argsort(idx)
    return idx[o], np.asarray(ys, float)[o], np.asarray(cs, float)[o]

def person_folds():
    rng = np.random.default_rng(C.SEED)
    f = np.arange(C.NTOT) % NFOLD
    rng.shuffle(f)
    return f

def foldstd(a, fold):
    out = np.zeros_like(a, dtype=np.float64)
    for f in range(NFOLD):
        tr = fold != f; te = fold == f
        m = a[tr].mean(); s = a[tr].std()
        out[te] = (a[te] - m) / (s if s > 1e-9 else 1.0)
    return out

def bh(p):
    p = np.asarray(p, float); m = np.isfinite(p)
    q = np.full(len(p), np.nan); pp = p[m]
    o = np.argsort(pp); n = len(pp)
    r = np.empty(n); r[o] = np.arange(1, n + 1)
    qs = np.minimum.accumulate((pp * n / r)[o][::-1])[::-1]
    out = np.empty(n); out[o] = qs
    q[m] = np.minimum(out, 1.0)
    return q

def main():
    idx, y, Xc = load_pheno()
    fold = person_folds()[idx]
    N = len(y)
    D = np.column_stack([np.ones(N), Xc])
    Q, _ = np.linalg.qr(D)
    yp = y - Q @ (Q.T @ y); yy = float(yp @ yp)
    dfree = N - D.shape[1] - 1
    print('N=%d, 열 %d개' % (N, len(COLS)), flush=True)

    def zof(S):
        Sp = S - Q @ (Q.T @ S)
        ss = float(Sp @ Sp)
        if ss <= 1e-8:
            return np.nan
        num = float(Sp @ yp)
        rss = yy - num * num / ss
        se = np.sqrt(max(rss, 0) / dfree / ss)
        return (num / ss) / se if se > 0 else np.nan

    Zs = {n: [] for n, _ in COLS}
    nd = 0
    for f in sorted(glob.glob(C.BURD + '/burden_chr*.npy')):
        ch = os.path.basename(f).replace('burden_chr', '').replace('.npy', '')
        dif = C.BURD + '/domidx_chr%s.npy' % ch
        if not os.path.exists(dif):
            continue
        dom = np.load(dif); A = np.load(f, mmap_mode='r')
        for i in range(len(dom)):
            blk = np.asarray(A[i][:, idx], dtype=np.float64)
            for nm, c in COLS:
                a = blk[c]
                Zs[nm].append(zof(foldstd(a, fold)) if a.std() > 1e-9
                              else np.nan)
            nd += 1
        del A
        print('  chr%-3s 누적 %d' % (ch, nd), flush=True)

    res = []
    zref = np.array(Zs['RAW_ALL'])
    for nm, _ in COLS:
        z = np.array(Zs[nm])
        q = bh(2.0 * norm.sf(np.abs(z)))
        m = np.isfinite(z) & np.isfinite(zref)
        res.append(dict(col=nm, n_dom=int(np.isfinite(z).sum()),
                        n_sig=int(np.nansum(q < 0.05)),
                        median_abs_z=float(np.nanmedian(np.abs(z))),
                        max_abs_z=float(np.nanmax(np.abs(z))),
                        mean_delta_vs_raw=float(
                            np.nanmean(np.abs(z[m]) - np.abs(zref[m]))),
                        corr_z_with_raw=float(
                            np.corrcoef(z[m], zref[m])[0, 1])))
    print()
    h = '%-12s %6s %7s %11s %10s %14s %12s'
    print(h % ('열', 'n', '유의', '|z|중앙', '|z|최대', 'Δ|z| vs RAW', 'z상관(RAW)'))
    for r in res:
        print(h % (r['col'], r['n_dom'], r['n_sig'],
                   '%.3f' % r['median_abs_z'], '%.2f' % r['max_abs_z'],
                   '%+.4f' % r['mean_delta_vs_raw'],
                   '%.3f' % r['corr_z_with_raw']))
    out = C.ATT + '/compare_std.json'
    json.dump(dict(n_domains=nd, trait='tchl', domstep=C.DOMSTEP,
                   result=res), open(out, 'w'), indent=1, ensure_ascii=False)
    print('\n->', out)

if __name__ == '__main__':
    main()
