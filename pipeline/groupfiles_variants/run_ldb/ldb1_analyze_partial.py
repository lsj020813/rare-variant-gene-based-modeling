#!/usr/bin/env python3
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
import numpy as np, json, sys, gzip, io
D = _configured('${PROJECT_ROOT}/work/run_ldb/out/partial')
rng = np.random.default_rng(20260914)
B = []
with open(f'{D}/blocks19.blocks.det') as fh:
    fh.readline()
    for line in fh:
        q = line.split()
        B.append((int(q[1]), int(q[2])))
B = np.array(sorted(B))
nb = len(B)
lens = B[:, 1] - B[:, 0]
chr_lo, chr_hi = (B[0, 0], B[-1, 1])

def block_id(pos, starts, ends):
    i = np.searchsorted(starts, pos, side='right') - 1
    ok = (i >= 0) & (pos <= ends[np.clip(i, 0, len(ends) - 1)])
    return np.where(ok, i, -1)
posA = []
posB = []
r2 = []
opener = gzip.open if f'{D}/r2_19.vcor'.endswith('.gz') else open
with open(f'{D}/r2_19.vcor') as fh:
    hdr = fh.readline().lstrip('#').split()
    ci = {h: i for i, h in enumerate(hdr)}
    ia, ib, ir = (ci['POS_A'], ci['POS_B'], [k for k in ci if 'R2' in k.upper()][0])
    ir = ci[ir]
    for line in fh:
        q = line.split()
        posA.append(int(q[ia]))
        posB.append(int(q[ib]))
        r2.append(float(q[ir]))
posA = np.array(posA)
posB = np.array(posB)
r2 = np.array(r2)
dist = np.abs(posA - posB)
bins = [(1000, 10000), (10000, 50000), (50000, 200000)]

def ratios(starts, ends):
    a = block_id(posA, starts, ends)
    b = block_id(posB, starts, ends)
    both = (a >= 0) & (b >= 0)
    same = both & (a == b)
    diff = both & (a != b)
    out = []
    for lo, hi in bins:
        m = (dist >= lo) & (dist < hi)
        ws = r2[m & same]
        wd = r2[m & diff]
        out.append(dict(n_same=int(len(ws)), n_diff=int(len(wd)), med_same=float(np.median(ws)) if len(ws) else None, med_diff=float(np.median(wd)) if len(wd) else None, ratio=float(np.median(ws) / np.median(wd)) if len(ws) and len(wd) and (np.median(wd) > 0) else None, frac_same_r2_gt0p1=float((ws > 0.1).mean()) if len(ws) else None, frac_diff_r2_gt0p1=float((wd > 0.1).mean()) if len(wd) else None))
    return out
obs = ratios(B[:, 0], B[:, 1])
null = [[] for _ in bins]
NP = 200
for it in range(NP):
    L = rng.permutation(lens)
    gaps = rng.dirichlet(np.ones(nb + 1)) * max(1, chr_hi - chr_lo - L.sum())
    s = chr_lo + np.cumsum(gaps[:-1]) + np.concatenate([[0], np.cumsum(L[:-1])])
    s = s.astype(int)
    e = s + L
    rr = ratios(s, e)
    for k in range(len(bins)):
        null[k].append(rr[k]['ratio'] if rr[k]['ratio'] is not None else np.nan)
res = dict(n_blocks=nb, block_len_kb_median=float(np.median(lens)) / 1000.0, block_cover_frac=float(lens.sum() / (chr_hi - chr_lo)), n_pairs=int(len(r2)), n_index=int(len(set(posA.tolist()))), bins=[dict(kb=f'{lo // 1000}-{hi // 1000}', **obs[k], null_ratio_median=float(np.nanmedian(null[k])), null_ratio_q95=float(np.nanquantile(null[k], 0.95)), pass_1p5=obs[k]['ratio'] is not None and obs[k]['ratio'] >= 1.5, pass_null=obs[k]['ratio'] is not None and obs[k]['ratio'] > float(np.nanquantile(null[k], 0.95))) for k, (lo, hi) in enumerate(bins)])
res['verdict_prereg'] = sum((1 for b in res['bins'] if b['pass_1p5'] and b['pass_null'])) >= 2
json.dump(res, open(f'{D}/ldb1_result.json', 'w'), indent=1)
print(json.dumps(res))
