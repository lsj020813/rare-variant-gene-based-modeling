#!/usr/bin/env python
import os as _config_os
import re as _config_re

def _required(name):
    value = _config_os.environ.get(name)
    if value is None or not value.strip():
        raise RuntimeError(f'Set {name} before running this script')
    return value

def _required_int(name):
    value = int(_required(name))
    if value <= 0:
        raise ValueError(f'{name} must be a positive integer')
    return value

def _configured(value):
    defaults = {'BCFTOOLS': 'bcftools', 'PLINK2': 'plink2', 'PYTHON': 'python3', 'TABIX': 'tabix', 'BGZIP': 'bgzip', 'SAMTOOLS': 'samtools', 'BEDTOOLS': 'bedtools', 'CROSSMAP': 'CrossMap', 'UDOCKER_BIN': 'udocker'}
    def resolve(match):
        name = match.group(1) or match.group(2)
        if name in _config_os.environ:
            return _required(name)
        return defaults.get(name) or _required(name)
    return _config_re.sub(r"\$\{([A-Z][A-Z0-9_]*)\}|\$([A-Z][A-Z0-9_]*)", resolve, value)

def _source_path(base, value):
    return value if _config_os.path.isabs(value) else _config_os.path.join(base, value)
import sys, json, time, numpy as np
from cyvcf2 import VCF
V = _configured('${PROJECT_ROOT}/work/ref/band_vcf/chr22.band.vcf.gz')
G = _configured('${PROJECT_ROOT}/work/ref/groupfiles_bwg/chr22.B_3kb_re2g.txt')
OUT = sys.argv[1]
t0 = time.time()
wins = {}
for line in open(G):
    p = line.split()
    if len(p) > 2 and p[1] == 'var':
        wins[p[0]] = p[2:]
vcf = VCF(V, gts012=True, lazy=True)

def fetch(vid):
    ch, pos, ref, alt = vid.split(':')
    for v in vcf(f'{ch}:{pos}-{pos}'):
        if v.REF == ref and v.ALT and (v.ALT[0] == alt):
            ds = np.asarray(v.format('DS'), dtype=np.float32).ravel()
            return (ds, float(v.INFO.get('R2', np.nan)), float(v.INFO.get('MAF', np.nan)), int(pos))
    return None
bins = np.array([0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 0.9, 0.95, 0.99, 1.0001])
strata = {}

def add(key, r2vals):
    h = np.histogram(r2vals, bins=bins)[0]
    strata[key] = strata.get(key, np.zeros(len(bins) - 1, dtype=np.int64)) + h

def meff(C):
    if C.shape[0] < 2:
        return float(C.shape[0])
    lam = np.linalg.eigvalsh(C)
    lam = np.clip(lam, 0, None)
    return float(lam.sum() ** 2 / (lam ** 2).sum())
perwin = []
nmiss = 0
nvar = 0
for gi, (g, vids) in enumerate(wins.items()):
    rows = []
    R = []
    MAF = []
    POS = []
    for vid in vids:
        f = fetch(vid)
        if f is None:
            nmiss += 1
            continue
        ds, r2, maf, pos = f
        rows.append(ds)
        R.append(r2)
        MAF.append(maf)
        POS.append(pos)
    M = len(rows)
    nvar += M
    if M < 2:
        perwin.append(dict(gene=g, M=M))
        continue
    X = np.vstack(rows)
    X = X - X.mean(1, keepdims=True)
    sd = X.std(1)
    sd[sd == 0] = np.nan
    C = X @ X.T / X.shape[1] / np.outer(sd, sd)
    C = np.nan_to_num(C)
    np.fill_diagonal(C, 1.0)
    R = np.array(R)
    POS = np.array(POS)
    iu = np.triu_indices(M, 1)
    r2p = C[iu] ** 2
    ra = R[iu[0]]
    rb = R[iu[1]]
    d = np.abs(POS[iu[0]] - POS[iu[1]])
    for thr in (0.8, 0.9):
        both = (ra >= thr) & (rb >= thr)
        one = (ra >= thr) ^ (rb >= thr)
        none = (ra < thr) & (rb < thr)
        for nm, mask in (('both', both), ('one', one), ('neither', none)):
            add(f'R2>={thr}:{nm}:all', r2p[mask])
            add(f'R2>={thr}:{nm}:d<3kb', r2p[mask & (d < 3000)])
            add(f'R2>={thr}:{nm}:d>=3kb', r2p[mask & (d >= 3000)])
    rec = dict(gene=g, M=M, meff_all=meff(C), maxr2_all=float(r2p.max()) if r2p.size else None, medR2=float(np.median(R)))
    for thr in (0.8, 0.9):
        k = R >= thr
        rec[f'n_R2ge{thr}'] = int(k.sum())
        rec[f'meff_R2ge{thr}'] = meff(C[np.ix_(k, k)]) if k.sum() >= 2 else float(k.sum())
        rec[f'maxr2_R2ge{thr}'] = float((C[np.ix_(k, k)][np.triu_indices(k.sum(), 1)] ** 2).max()) if k.sum() >= 2 else None
    perwin.append(rec)
    if gi % 40 == 0:
        print(f'progress {gi}/{len(wins)} nvar={nvar} t={time.time() - t0:.0f}s', flush=True)
out = dict(n_windows=len(wins), n_var_loaded=nvar, n_var_missing=nmiss, bins=bins.tolist(), strata={k: v.tolist() for k, v in strata.items()}, perwin=perwin, elapsed_s=time.time() - t0)
json.dump(out, open(OUT, 'w'))
print('DONE', json.dumps({k: int(v.sum()) for k, v in strata.items() if k.endswith(':all')}))
