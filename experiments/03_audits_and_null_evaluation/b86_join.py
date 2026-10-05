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
import os, json, glob, numpy as np
from scipy.stats import norm
R15 = _configured('${PROJECT_ROOT}/work/ref15')
OUT = R15 + '/saige_single_typed'
W = _configured('${PROJECT_ROOT}/work')
traits = ['tchl', 'htn', 'dm', 'lip']

def load_groups(N):
    g = {}
    for line in open(f'{R15}/groupfiles_bwg/chr{N}.B_3kb_re2g.txt'):
        q = line.split()
        if len(q) > 2 and q[1] == 'var':
            g[q[0]] = q[2:]
    return g
P = []
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
        for gene, vids in groups.items():
            zs = [abs(z[v]) for v in vids if v in z]
            ps = [pv[v] for v in vids if v in pv]
            if len(zs) < 2:
                continue
            pos = [int(v.split(':')[1]) for v in vids if v in z]
            P.append((T, N, gene, len(zs), min(ps), max(zs), min(pos), max(pos)))
n = len(P)
pmin = np.array([p[4] for p in P])
M = np.array([p[3] for p in P])
psid = 1 - (1 - pmin) ** M
order = np.argsort(psid)
q = 0.2
kmax = max([k for k, i in enumerate(order, 1) if psid[i] <= q * k / n] or [0])
cand = [P[i] for i in order[:kmax]]
assert len(cand) == 86, len(cand)
regions = {}
for f in glob.glob(f'{W}/run_ourfm/regions/*.regions.tsv'):
    T, ch = os.path.basename(f).split('.')[0:2]
    N = int(ch[3:])
    with open(f) as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        for line in fh:
            q = dict(zip(hdr, line.rstrip('\n').split('\t')))
            s = int(float(q.get('start', q.get('region_start', 0))))
            e = int(float(q.get('end', q.get('region_end', 0))))
            regions.setdefault((T, N), []).append((s, e))
lead = {}
for T in traits:
    for N in range(1, 23):
        f = f'{W}/ref/gwas05/{T}.chr{N}.txt'
        if not os.path.exists(f):
            continue
        with open(f) as fh:
            hdr = fh.readline().rstrip('\n').split('\t')
            ci = {h: i for i, h in enumerate(hdr)}
            L = []
            for line in fh:
                q = line.split('\t')
                try:
                    p = float(q[ci['p.value']])
                except:
                    continue
                if p < 5e-08:
                    L.append(int(q[ci['POS']]))
            lead[T, N] = np.array(sorted(L))
rows = []
for T, N, g, Mw, pm, mz, lo, hi in cand:
    in_reg = any((s - 500000 <= hi and lo <= e + 500000 for s, e in regions.get((T, N), [])))
    L = lead.get((T, N), np.array([]))
    near = bool(L.size) and bool(((L >= lo - 500000) & (L <= hi + 500000)).any())
    any_trait = any((bool(lead.get((T2, N), np.array([])).size) and bool(((lead[T2, N] >= lo - 500000) & (lead[T2, N] <= hi + 500000)).any()) for T2 in traits))
    rows.append(dict(trait=T, chr=N, gene=g, M=Mw, pmin=pm, maxz=round(mz, 2), lo=lo, hi=hi, in_region=in_reg, near_lead_same_trait=near, near_lead_any_trait=any_trait))
summ = dict(n_cand=len(rows), in_region=sum((r['in_region'] for r in rows)), near_lead_same=sum((r['near_lead_same_trait'] for r in rows)), near_lead_any=sum((r['near_lead_any_trait'] for r in rows)), novel_strict=sum((not r['near_lead_any_trait'] for r in rows)), novel_by_trait={T: sum((1 for r in rows if r['trait'] == T and (not r['near_lead_any_trait']))) for T in traits}, novel_list=[(r['trait'], r['chr'], r['gene'], r['M'], r['pmin'], r['maxz']) for r in rows if not r['near_lead_any_trait']], regions_loaded=sum((len(v) for v in regions.values())), lead_files=len(lead), n_regions_files=len(regions))
json.dump(dict(summary=summ, rows=rows), open(f'{W}/run_meth/out/gate2_B86_vs_gwas.json', 'w'), indent=1)
print(json.dumps(summ, default=str))
