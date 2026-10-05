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
import json, os, numpy as np
W = _configured('${PROJECT_ROOT}/work')
traits = ['tchl', 'htn', 'dm', 'lip']
K = json.load(open(f'{W}/run_meth/known_loci_grch37.json'))
d = json.load(open(f'{W}/run_meth/out/gate2_B86_vs_gwas.json'))
cand = d['rows']

def load_groups(N):
    g = {}
    for line in open(f'{W}/ref15/groupfiles_bwg/chr{N}.B_3kb_re2g.txt'):
        q = line.split()
        if len(q) > 2 and q[1] == 'var':
            g[q[0]] = q[2:]
    return g
lead = {}
for T in traits:
    for N in range(1, 23):
        f = f'{W}/ref/gwas05/{T}.chr{N}.txt'
        L = []
        with open(f) as fh:
            hdr = fh.readline().rstrip('\n').split('\t')
            ci = {h: i for i, h in enumerate(hdr)}
            for line in fh:
                q = line.split('\t')
                try:
                    p = float(q[ci['p.value']])
                except:
                    continue
                if p < 5e-08:
                    L.append(int(q[ci['POS']]))
        lead[T, N] = np.array(sorted(L))
ext = {}
for N in range(1, 23):
    f = f'{W}/ref/gwas05/tchl.chr{N}.txt'
    with open(f) as fh:
        fh.readline()
        first = int(fh.readline().split('\t')[1])
    last = int(os.popen(f'tail -1 {f} | cut -f2').read().strip())
    ext[N] = (first, last)
FDR = {T: {} for T in traits}
for r in cand:
    FDR[r['trait']].setdefault(r['chr'], []).append((r['lo'], r['hi'], r['M'], r['pmin']))
FW = 0.05 / 27024

def overlaps(T, N, lo, hi, fwer=False):
    for a, b, M, pm in FDR[T].get(N, []):
        if a <= hi and b >= lo:
            if not fwer or 1 - (1 - pm) ** M < FW:
                return True
    return False
res = {}
rng = np.random.default_rng(20260914)
for T in traits:
    loci = [k for k in K if k['trait'] == T]
    a = b = 0
    kF = kW = 0
    bl = []
    hits = []
    for k in loci:
        N = k['chr']
        lo, hi = (k['lo'] - 500000, k['hi'] + 500000)
        L = lead[T, N]
        ours = bool(L.size) and bool(((L >= lo) & (L <= hi)).any())
        if ours:
            a += 1
            continue
        b += 1
        bl.append((N, lo, hi))
        if overlaps(T, N, lo, hi):
            kF += 1
            hits.append((N, k['lo'], k['hi'], k['n']))
        if overlaps(T, N, lo, hi, fwer=True):
            kW += 1
    null = np.zeros(1000, int)
    for it in range(1000):
        cnt = 0
        for N, lo, hi in bl:
            w = hi - lo
            s = rng.integers(ext[N][0], max(ext[N][0] + 1, ext[N][1] - w))
            if overlaps(T, N, s, s + w):
                cnt += 1
        null[it] = cnt
    res[T] = dict(known_loci=len(loci), a_ours_ge=a, b_missed=b, frac_a=round(a / max(1, len(loci)), 3), k_FDR=kF, k_FWER=kW, null_mean=round(float(null.mean()), 2), null_q95=int(np.quantile(null, 0.95)), emp_p=round(float((null >= kF).mean()), 3), hits=hits[:12])
tot = dict(k_FDR=sum((r['k_FDR'] for r in res.values())), k_FWER=sum((r['k_FWER'] for r in res.values())), b=sum((r['b_missed'] for r in res.values())), null_mean=round(sum((r['null_mean'] for r in res.values())), 2))
json.dump(dict(per_trait=res, total=tot, rule='prereg v151'), open(f'{W}/run_meth/out/complement_test.json', 'w'), indent=1)
print(json.dumps(dict(total=tot, per_trait={T: {k: v for k, v in r.items() if k != 'hits'} for T, r in res.items()})))
print('HITS', json.dumps({T: r['hits'] for T, r in res.items()}))
