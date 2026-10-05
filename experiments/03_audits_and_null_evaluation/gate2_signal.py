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
import os, sys, json, math, glob, numpy as np
from scipy.stats import norm
R15 = _configured('${PROJECT_ROOT}/work/ref15')
OUT = R15 + '/saige_single_typed'
traits = ['tchl', 'htn', 'dm', 'lip']
p3 = 2 * norm.sf(3.0)

def load_groups(N):
    g = {}
    for line in open(f'{R15}/groupfiles_bwg/chr{N}.B_3kb_re2g.txt'):
        q = line.split()
        if len(q) > 2 and q[1] == 'var':
            g[q[0]] = q[2:]
    return g
res = {'arms_done': [], 'pairs': []}
for N in range(1, 23):
    groups = None
    for T in traits:
        f = f'{OUT}/{T}.chr{N}'
        if not os.path.exists(f + '.done'):
            continue
        if groups is None:
            groups = load_groups(N)
        z = {}
        pv = {}
        with open(f) as fh:
            hdr = fh.readline().rstrip('\n').split('\t')
            ci = {h: i for i, h in enumerate(hdr)}
            for line in fh:
                q = line.rstrip('\n').split('\t')
                vid = f"{q[ci['CHR']]}:{q[ci['POS']]}:{q[ci['Allele1']]}:{q[ci['Allele2']]}"
                try:
                    b = float(q[ci['BETA']])
                    s = float(q[ci['SE']])
                    p = float(q[ci['p.value']])
                except:
                    continue
                if s > 0:
                    z[vid] = b / s
                    pv[vid] = p
        res['arms_done'].append(f'{T}.chr{N}')
        for gene, vids in groups.items():
            zs = [abs(z[v]) for v in vids if v in z]
            ps = [pv[v] for v in vids if v in pv]
            if len(zs) < 2:
                continue
            res['pairs'].append((T, N, gene, len(zs), min(ps), max(zs)))
P = res['pairs']
n = len(P)
if n == 0:
    print('NO_PAIRS_YET')
    sys.exit(0)
M = np.array([p[3] for p in P])
pmin = np.array([p[4] for p in P])
maxz = np.array([p[5] for p in P])
N0 = np.sum(1 - (1 - p3) ** M)
N0_sd = math.sqrt(np.sum((1 - (1 - p3) ** M) * (1 - p3) ** M))
Nobs = int((maxz >= 3).sum())
excess = Nobs - N0
psid = 1 - (1 - pmin) ** M
order = np.argsort(psid)
q = 0.2
thr = [psid[i] for k, i in enumerate(order, 1) if psid[i] <= q * k / n]
nB = len(thr)
n45 = int((maxz >= 4.5).sum())
n545 = int((maxz >= 5.45).sum())
out = dict(arms_done=len(res['arms_done']), pairs=n, windows_typed_ge2=len(set(((p[1], p[2]) for p in P))), A_Nobs_z3=Nobs, A_N0=round(N0, 1), A_N0_sd=round(N0_sd, 1), A_excess=round(excess, 1), A_excess_z=round(excess / N0_sd, 2) if N0_sd > 0 else None, B_fdr20_pairs=nB, ref_maxz_ge4_5=n45, ref_maxz_ge5_45=n545, by_trait={T: dict(pairs=int(sum((1 for p in P if p[0] == T))), z3=int(sum((1 for p in P if p[0] == T and p[5] >= 3))), z45=int(sum((1 for p in P if p[0] == T and p[5] >= 4.5)))) for T in traits})
print(json.dumps(out))
json.dump(out, open(sys.argv[1] if len(sys.argv) > 1 else '/dev/null', 'w'))
