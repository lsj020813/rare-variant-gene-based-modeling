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
NPERM = int(os.environ.get('NPERM', '20'))
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
        m = a[tr].mean(); s = a[tr].std()
        out[te] = (a[te] - m) / (s if s > 1e-9 else 1.0)
    return out

def main():
    Z = np.load(C.LABELS, allow_pickle=True)
    names = [str(x) for x in Z['names']]
    kinds = [str(x) for x in Z['kinds']]
    P = len(names)

    idx, y, Xc = load_pheno()
    fold = person_folds()[idx]
    N = len(y)
    D = np.column_stack([np.ones(N), Xc])
    Q, _ = np.linalg.qr(D)

    yp0 = y - Q @ (Q.T @ y)
    rng = np.random.default_rng(C.SEED + 5150)
    YP = np.empty((N, NPERM + 1))
    YP[:, 0] = yp0
    for k in range(1, NPERM + 1):
        YP[:, k] = rng.permutation(yp0)
    YY = (YP * YP).sum(0)
    dfree = N - D.shape[1] - 1
    print('N=%d, 순열 %d회, w 격자 %d점' % (N, NPERM, len(WGRID)), flush=True)

    def zvec(S):
        Sp = S - Q @ (Q.T @ S)
        ss = float(Sp @ Sp)
        if ss <= 1e-8:
            return None
        num = Sp @ YP
        rss = YY - num * num / ss
        se = np.sqrt(np.maximum(rss, 0) / dfree / ss)
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.where(se > 0, (num / ss) / se, np.nan)

    G = {n: [] for n in names}
    CO = {n: [] for n in names}
    ZF = []
    nd = 0
    for f in sorted(glob.glob(C.BURD + '/burden_chr*.npy')):
        ch = os.path.basename(f).replace('burden_chr', '').replace('.npy', '')
        dif = C.BURD + '/domidx_chr%s.npy' % ch
        if not os.path.exists(dif):
            continue
        dom = np.load(dif)
        A = np.load(f, mmap_mode='r')
        for i in range(len(dom)):
            blk = np.asarray(A[i][:, idx], dtype=np.float64)
            zf = zvec(foldstd(blk[C.C_TOT_ALL], fold))
            if zf is None or not np.isfinite(zf[0]):
                continue
            ZF.append(abs(zf[0])); nd += 1
            for p, nm in enumerate(names):
                a = blk[C.C_P0 + 2 * p]; b = blk[C.C_P0 + 2 * p + 1]
                if a.std() < 1e-9 or b.std() < 1e-9:
                    G[nm].append(np.full(NPERM + 1, np.nan))
                    CO[nm].append(np.nan)
                    continue
                CO[nm].append(float(np.corrcoef(a, b)[0, 1]))
                za = foldstd(a, fold); zb = foldstd(b, fold)
                Zs = []
                for w in WGRID:
                    v = zvec(w * za + (1 - w) * zb)
                    if v is not None:
                        Zs.append(np.abs(v))
                if not Zs:
                    G[nm].append(np.full(NPERM + 1, np.nan)); continue
                G[nm].append(np.nanmax(np.vstack(Zs), 0) - np.abs(zf))
        del A
        print('  chr%-3s 누적 %d' % (ch, nd), flush=True)

    ZF = np.array(ZF)
    res = []
    for p, nm in enumerate(names):
        M = np.vstack(G[nm])
        med = np.nanmedian(M, 0)
        obs, nul = med[0], med[1:]
        pval = (1 + int((nul >= obs).sum())) / (NPERM + 1)
        res.append(dict(
            partition=nm, kind=kinds[p], n_dom=int(np.isfinite(M[:, 0]).sum()),
            corr_median=float(np.nanmedian(CO[nm])),
            gain_obs=float(obs), gain_null_mean=float(np.nanmean(nul)),
            gain_null_lo=float(np.nanmin(nul)), gain_null_hi=float(np.nanmax(nul)),
            excess=float(obs - np.nanmean(nul)), p=pval))
    res.sort(key=lambda r: -r['excess'])
    print()
    h = '%-15s %-17s %8s %9s %10s %16s %7s'
    print(h % ('분할', '축', '상관중앙', '관측이득', '귀무평균', '귀무범위', '순열p'))
    for r in res:
        print(h % (r['partition'], r['kind'], '%.3f' % r['corr_median'],
                   '%.4f' % r['gain_obs'], '%.4f' % r['gain_null_mean'],
                   '[%.4f,%.4f]' % (r['gain_null_lo'], r['gain_null_hi']),
                   '%.3f' % r['p']))
    out = C.ATT + '/ceiling_cand2.json'
    json.dump(dict(n_domains=nd, w_grid=list(WGRID), nperm=NPERM,
                   trait='tchl', domstep=C.DOMSTEP,
                   p_floor=1.0 / (NPERM + 1), result=res),
              open(out, 'w'), indent=1, ensure_ascii=False)
    print('\n->', out)

if __name__ == '__main__':
    main()
