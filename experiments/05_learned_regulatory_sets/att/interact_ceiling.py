import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import json
import os
import sys
import time
import numpy as np
from scipy.linalg import cho_factor, cho_solve

sys.path.insert(0, _config_path('${PROJECT_ROOT}/work/fset/att'))
import att_common as C
import att_common_oracle as CO

CACHE = C.FSET + '/att/hqcache'
OUT = C.FSET + '/att/interact'
MAXV = int(os.environ.get('MAXV', '50'))
MINV = int(os.environ.get('MINV', '8'))
ALPHAS = CO.ALPHAS
NFOLD = CO.NFOLD

def oof_z(X, PH):
    Y = PH['Y']
    bnd = PH['bounds']
    P = X.shape[1]
    I = np.eye(P)
    Gf = X.T @ X
    XtYf = X.T @ Y
    best = None
    for a in ALPHAS:
        pred = np.zeros_like(Y)
        ok = True
        for f in range(NFOLD):
            a0, a1 = bnd[f], bnd[f + 1]
            Xte = X[a0:a1]
            G = Gf - Xte.T @ Xte
            b = XtYf - Xte.T @ Y[a0:a1]
            ntr = X.shape[0] - (a1 - a0)
            try:
                c = cho_factor(G + (a * ntr) * I, lower=True)
                pred[a0:a1] = Xte @ cho_solve(c, b)
            except Exception:
                ok = False
                break
        if not ok:
            continue
        zs = np.full(Y.shape[1], np.nan)
        for k in range(Y.shape[1]):
            v = CO._z(pred[:, k], PH)
            if v is not None:
                zs[k] = abs(v[k])
        best = zs if best is None else np.fmax(best, zs)
    return best

def main():
    ch = sys.argv[1]
    os.makedirs(OUT, exist_ok=True)
    idxf = CACHE + '/hqidx_chr%s.npy' % ch
    if not os.path.exists(idxf):
        print('캐시 없음 chr%s' % ch)
        return
    vidx = np.load(idxf)
    got = np.load(CACHE + '/hqgot_chr%s.npy' % ch)
    A = np.load(CACHE + '/hq_chr%s.npy' % ch, mmap_mode='r')
    d = C.master()
    di = np.asarray(d['dom_idx'])
    PH = CO.pheno()
    pidx, fold = PH['idx'], PH['fold']
    rng = np.random.default_rng(C.SEED + 31415)

    doms = sorted(set(di[vidx[got]].tolist()))
    print('chr%s 도메인 %d개' % (ch, len(doms)), flush=True)
    rep = []
    t0 = time.time()
    for g in doms:
        rows = np.where((di[vidx] == g) & got)[0]
        if len(rows) < MINV:
            continue
        if len(rows) > MAXV:
            rows = np.sort(rng.choice(rows, MAXV, replace=False))
        X0 = np.asarray(A[rows][:, pidx], dtype=np.float64).T
        sd = X0.std(0)
        keep = sd > 1e-9
        X0 = X0[:, keep]
        M = X0.shape[1]
        if M < MINV:
            continue
        Xs = np.zeros_like(X0)
        for f in range(NFOLD):
            tr = fold != f
            te = fold == f
            m = X0[tr].mean(0)
            s = X0[tr].std(0)
            s[s < 1e-9] = 1.0
            Xs[te] = (X0[te] - m) / s
        del X0
        iu, ju = np.triu_indices(M, k=1)
        XI = Xs[:, iu] * Xs[:, ju]
        XI = XI - XI.mean(0)
        sdI = XI.std(0)
        sdI[sdI < 1e-9] = 1.0
        XI /= sdI
        zm = oof_z(Xs, PH)
        zi = oof_z(np.hstack([Xs, XI]), PH)
        r = dict(domain_idx=int(g), chrom=str(ch), n_var=int(M),
                 n_pairs=int(len(iu)))
        for nm, z in (('MAIN', zm), ('MAIN_INT', zi)):
            if z is None:
                continue
            r[nm] = dict(obs=float(z[0]), null_mean=float(np.nanmean(z[1:])),
                         null_max=float(np.nanmax(z[1:])),
                         excess=float(z[0] - np.nanmean(z[1:])),
                         p=float((1 + int((z[1:] >= z[0]).sum())) / len(z)))
        rep.append(r)
        del Xs, XI
    json.dump(dict(chrom=str(ch), n_domains=len(rep), maxv=MAXV,
                   wall_s=round(time.time() - t0, 1), per_domain=rep),
              open(OUT + '/int_chr%s.json' % ch, 'w'))
    open(OUT + '/int_chr%s.done' % ch, 'w').write(
        '%d domains %.0f s\n' % (len(rep), time.time() - t0))
    print('chr%s 완료 %d 도메인 %.0f초' % (ch, len(rep), time.time() - t0),
          flush=True)

if __name__ == '__main__':
    main()
