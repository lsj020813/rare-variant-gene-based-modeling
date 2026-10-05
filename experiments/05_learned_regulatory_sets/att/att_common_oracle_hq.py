import os as _cfg_os
import re as _cfg_re
import math as _cfg_math

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)

def _config_number(name, cast=int, positive=False):
    value = cast(_cfg_os.environ[name])
    if not _cfg_math.isfinite(value) or (positive and value <= 0):
        raise ValueError(name + " must be finite" + (" and positive" if positive else ""))
    return value

import os as _os
N_SAMPLES = _config_number("N_SAMPLES", int, True)
import csv
import os
import numpy as np
from scipy.linalg import cho_factor, cho_solve

FSET = _config_path('${PROJECT_ROOT}/work/fset')
ATT = FSET + '/att/oracle_hq'
PRIV = ATT + '/private'
BURD = PRIV + '/burden'
LOGS = ATT + '/logs'
MASTER = FSET + '/f3f6/cache/master.npz'
VCFD = _config_path('${PROJECT_ROOT}/work/ref/orig_index')
PHENO = _config_path('${PROJECT_ROOT}/work/ref/pheno_v3/tchl_v3.tsv')
LABELS = FSET + '/att/cand_labels.npz'

SEED = 20260923
NTOT = N_SAMPLES
NDOM = 1000
NFOLD = 5
NPERM = int(os.environ.get('NPERM', '10'))
ALPHAS = [0.01, 0.1, 1.0, 10.0]
COVARS = ['age', 'sex_male', 'CT', 'NC', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5']
STRATA = [
    ('ALL',       0.0,  1.01, 0.0, 1.01),
    ('HQ_ALL',    0.0,  1.01, 0.9, 1.01),
    ('HQ_RARE',   0.0,  0.01, 0.9, 1.01),
    ('LQ_RARE',   0.0,  0.01, 0.0, 0.5),
    ('HQ_COMMON', 0.05, 1.01, 0.9, 1.01),
]
MAF_SEC = 0.01
CARRIER_DS = 0.5
DOMSTEP = int(os.environ.get('DOMSTEP', '1'))

_PH = None

def master():
    return np.load(MASTER, allow_pickle=True)

def ensure_dirs():
    for p in (ATT, PRIV, BURD, LOGS):
        os.makedirs(p, exist_ok=True)

def pheno():
    global _PH
    if _PH is not None:
        return _PH
    from cyvcf2 import VCF
    v = VCF(VCFD + '/chr22.vcf.gz')
    pos = dict((s, i) for i, s in enumerate(v.samples))
    idx, ys, cs = [], [], []
    for r in csv.DictReader(open(PHENO), delimiter='\t'):
        i = pos.get(r['sample_id'])
        if i is None:
            continue
        try:
            yv = float(r['y']); cv = [float(r[c]) for c in COVARS]
        except (ValueError, TypeError, KeyError):
            continue
        if not np.isfinite(yv) or not np.isfinite(cv).all():
            continue
        idx.append(i); ys.append(yv); cs.append(cv)
    idx = np.asarray(idx, dtype=np.int64); o = np.argsort(idx)
    idx = idx[o]
    y = np.asarray(ys, float)[o]
    Xc = np.asarray(cs, float)[o]
    n = len(y)
    D = np.column_stack([np.ones(n), Xc])
    Q, _ = np.linalg.qr(D)
    yp = y - Q @ (Q.T @ y)
    rng = np.random.default_rng(SEED + 4242)
    Y = np.empty((n, NPERM + 1))
    Y[:, 0] = yp
    for k in range(1, NPERM + 1):
        Y[:, k] = rng.permutation(yp)
    Y = (Y - Y.mean(0)) / Y.std(0)
    fr = np.random.default_rng(SEED)
    f = np.arange(NTOT) % NFOLD
    fr.shuffle(f)
    fold = f[idx]
    order = np.argsort(fold, kind='stable')
    idx = idx[order]; Y = Y[order]; Q = Q[order]; fold = fold[order]
    b = [0]
    for k in range(NFOLD):
        b.append(int(np.searchsorted(fold, k, side='right')))
    _PH = dict(idx=idx, Y=Y, fold=fold, Q=Q, n=n,
               bounds=b, dfree=n - D.shape[1] - 1,
               YY=(Y * Y).sum(0))
    return _PH

def _z(S, PH):
    Q = PH['Q']
    Sp = S - Q @ (Q.T @ S)
    ss = float(Sp @ Sp)
    if ss <= 1e-8:
        return None
    num = Sp @ PH['Y']
    rss = PH['YY'] - num * num / ss
    se = np.sqrt(np.maximum(rss, 0) / PH['dfree'] / ss)
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(se > 0, (num / ss) / se, np.nan)

def oracle_domain(Dmat, got, mm, rr, g, ch, PH):
    idx, fold, Y = PH['idx'], PH['fold'], PH['Y']
    X0 = np.asarray(Dmat[idx], dtype=np.float64)
    out = dict(domain_idx=g, chrom=ch, n_var=int(Dmat.shape[1]))
    for nm, lo, hi, rlo, rhi in STRATA:
        sel = (got & (mm >= lo) & (mm < hi)
               & (rr >= rlo) & (rr < rhi) & np.isfinite(rr))
        sel = sel & (X0.std(0) > 1e-9)
        M = int(sel.sum())
        rec = dict(n=M)
        if M < 2:
            out[nm] = rec
            continue
        X = X0[:, sel]
        X = (X - X.mean(0)) / X.std(0)
        zf = _z(X.sum(1), PH)
        rec['flat_z'] = None if zf is None else float(abs(zf[0]))
        rec['flat_z_null'] = None if zf is None else float(
            np.nanmax(np.abs(zf[1:])))
        num = X.T @ Y
        ss = (X * X).sum(0)[:, None]
        rs = PH['YY'][None, :] - num * num / ss
        se = np.sqrt(np.maximum(rs, 0) / PH['dfree'] / ss)
        with np.errstate(divide='ignore', invalid='ignore'):
            zsing = np.where(se > 0, (num / ss) / se, np.nan)
        rec['single_z'] = float(np.nanmax(np.abs(zsing[:, 0])))
        rec['single_z_null_mean'] = float(np.nanmean(
            np.nanmax(np.abs(zsing[:, 1:]), 0)))
        Ps = {a: np.zeros_like(Y) for a in ALPHAS}
        okf = True
        I = np.eye(M)
        bnd = PH['bounds']
        Gfull = X.T @ X
        XtYfull = X.T @ Y
        for f in range(NFOLD):
            a0, a1 = bnd[f], bnd[f + 1]
            Xte = X[a0:a1]
            G = Gfull - Xte.T @ Xte
            XtY = XtYfull - Xte.T @ Y[a0:a1]
            ntr = X.shape[0] - (a1 - a0)
            for a in ALPHAS:
                try:
                    c = cho_factor(G + (a * ntr) * I, lower=True)
                    Ps[a][a0:a1] = Xte @ cho_solve(c, XtY)
                except Exception:
                    okf = False
            if not okf:
                break
        best = None
        if okf:
            for a in ALPHAS:
                zs = np.full(Y.shape[1], np.nan)
                for k in range(Y.shape[1]):
                    v = _z(Ps[a][:, k], PH)
                    if v is not None:
                        zs[k] = abs(v[k])
                best = zs if best is None else np.fmax(best, zs)
        if best is not None:
            rec['oof_z'] = float(best[0])
            rec['oof_z_null_mean'] = float(np.nanmean(best[1:]))
            rec['oof_z_null_max'] = float(np.nanmax(best[1:]))
            rec['oof_p'] = float((1 + int((best[1:] >= best[0]).sum()))
                                 / (len(best)))
        out[nm] = rec
    return out
