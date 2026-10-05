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

sys.path.insert(0, _config_path('${PROJECT_ROOT}/work/fset/att'))
import att_common_cand as C

WGRID = np.array([float(x) for x in os.environ.get(
    'WGRID', '-0.5,-0.25,0,0.25,0.5,0.75,1,1.25,1.5').split(',')])
NFOLD = 5

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
            yv = float(r['y'])
            cv = [float(r[c]) for c in C.COVARS]
        except (ValueError, TypeError, KeyError):
            continue
        if not np.isfinite(yv) or not np.isfinite(cv).all():
            continue
        idx.append(i); ys.append(yv); cs.append(cv)
    idx = np.asarray(idx, dtype=np.int64)
    o = np.argsort(idx)
    return (idx[o], np.asarray(ys, float)[o], np.asarray(cs, float)[o])

def person_folds():
    rng = np.random.default_rng(C.SEED)
    f = np.arange(C.NTOT) % NFOLD
    rng.shuffle(f)
    return f

def foldstd(a, fold):
    out = np.zeros_like(a, dtype=np.float64)
    for f in range(NFOLD):
        tr = fold != f
        te = fold == f
        m = a[tr].mean()
        s = a[tr].std()
        out[te] = (a[te] - m) / (s if s > 1e-9 else 1.0)
    return out

def main():
    Z = np.load(C.LABELS, allow_pickle=True)
    names = [str(x) for x in Z['names']]
    kinds = [str(x) for x in Z['kinds']]
    P = len(names)
    assert P == C.NPART

    idx, y, Xc = load_pheno()
    fold_all = person_folds()
    fold = fold_all[idx]
    D = np.column_stack([np.ones(len(y)), Xc])
    Q, _ = np.linalg.qr(D)
    yp = y - Q @ (Q.T @ y)
    yy = float(yp @ yp)
    dfree = len(y) - D.shape[1] - 1

    def zof(S):
        Sp = S - Q @ (Q.T @ S)
        ss = float(Sp @ Sp)
        if ss <= 1e-8:
            return np.nan
        num = float(Sp @ yp)
        beta = num / ss
        rss = yy - num * num / ss
        se = np.sqrt(rss / dfree / ss)
        return beta / se if se > 0 else np.nan

    gains = {n: [] for n in names}
    spans = {n: [] for n in names}
    corrs = {n: [] for n in names}
    zflats = []
    ndom = 0

    for f in sorted(glob.glob(C.BURD + '/burden_chr*.npy')):
        ch = os.path.basename(f).replace('burden_chr', '').replace('.npy', '')
        dif = C.BURD + '/domidx_chr%s.npy' % ch
        if not os.path.exists(dif):
            continue
        dom = np.load(dif)
        A = np.load(f, mmap_mode='r')
        for i in range(len(dom)):
            blk = np.asarray(A[i][:, idx], dtype=np.float64)
            zf = zof(foldstd(blk[C.C_TOT_ALL], fold))
            if not np.isfinite(zf):
                continue
            zflats.append(abs(zf))
            ndom += 1
            for p, nm in enumerate(names):
                a = blk[C.C_P0 + 2 * p]
                b = blk[C.C_P0 + 2 * p + 1]
                if a.std() < 1e-9 or b.std() < 1e-9:
                    gains[nm].append(np.nan)
                    spans[nm].append(np.nan)
                    corrs[nm].append(np.nan)
                    continue
                corrs[nm].append(float(np.corrcoef(a, b)[0, 1]))
                za = foldstd(a, fold)
                zb = foldstd(b, fold)
                zs = [abs(zof(w * za + (1 - w) * zb)) for w in WGRID]
                zs = [v for v in zs if np.isfinite(v)]
                if not zs:
                    gains[nm].append(np.nan); spans[nm].append(np.nan); continue
                gains[nm].append(max(zs) - abs(zf))
                spans[nm].append(max(zs) - min(zs))
        del A
        print('  chr%-3s 누적 도메인 %d' % (ch, ndom), flush=True)

    zflats = np.array(zflats)
    top = zflats >= np.percentile(zflats, 95)
    res = []
    for p, nm in enumerate(names):
        g = np.array(gains[nm], dtype=float)
        s = np.array(spans[nm], dtype=float)
        c = np.array(corrs[nm], dtype=float)
        m = np.isfinite(g)
        res.append(dict(
            partition=nm, kind=kinds[p], n_dom=int(m.sum()),
            corr_median=float(np.nanmedian(c)),
            gain_median=float(np.nanmedian(g)),
            gain_p90=float(np.nanpercentile(g, 90)),
            gain_top5_median=float(np.nanmedian(g[m & top[:len(g)]]))
            if (m & top[:len(g)]).any() else float('nan'),
            frac_gain_gt_0p5=float(np.nanmean(g > 0.5)),
            span_median=float(np.nanmedian(s))))

    res.sort(key=lambda r: -r['gain_median'])
    print()
    hdr = '%-16s %-17s %6s %8s %9s %9s %11s %9s'
    print(hdr % ('분할', '축', 'n', '상관중앙', '이득중앙', '이득p90',
                 '이득(상위5%)', '>0.5비율'))
    for r in res:
        print(hdr % (r['partition'], r['kind'], r['n_dom'],
                     '%.3f' % r['corr_median'], '%.4f' % r['gain_median'],
                     '%.4f' % r['gain_p90'], '%.4f' % r['gain_top5_median'],
                     '%.1f%%' % (100 * r['frac_gain_gt_0p5'])))
    out = C.ATT + '/ceiling_cand.json'
    json.dump(dict(n_domains=ndom, w_grid=list(WGRID), trait='tchl',
                   domstep=C.DOMSTEP, result=res),
              open(out, 'w'), indent=1, ensure_ascii=False)
    print('\n->', out)

if __name__ == '__main__':
    main()
