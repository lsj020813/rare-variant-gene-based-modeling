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
import glob, gzip, collections, random, statistics as st
R = _configured('${PROJECT_ROOT}/work/ref')

def read(files):
    out = {}
    for fp in files:
        with open(fp) as fh:
            h = fh.readline().rstrip('\n').split('\t')
            if 'Pvalue' not in h:
                continue
            gi = h.index('Region')
            pi = h.index('Pvalue')
            ni = h.index('Number_rare')
            for line in fh:
                fl = line.rstrip('\n').split('\t')
                try:
                    out[fl[gi]] = (float(fl[pi]), int(float(fl[ni])))
                except (ValueError, IndexError):
                    pass
    return out
B = read(sorted(glob.glob(f'{R}/saige_step2_pilot/tchl.chr19.B_3kb_re2g.B.part[0-9][0-9][0-9]')))
D = read([f'{R}/saige_step2_pilot/tchl.chr19.D_3kb_only'])
print(f'arm B genes: {len(B)}   arm D genes: {len(D)}')
import math
pairs = [(n, p) for p, n in B.values() if n > 0 and 0 < p <= 1]
bins = collections.defaultdict(list)
for n, p in pairs:
    b = 0 if n < 25 else 1 if n < 50 else 2 if n < 80 else 3 if n < 150 else 4
    bins[b].append(p)
lbl = {0: '<25', 1: '25-49', 2: '50-79', 3: '80-149', 4: '>=150'}
print('\n=== arm B: does p track variant count? ===')
for b in sorted(bins):
    v = bins[b]
    print(f'  n {lbl[b]:<8} genes={len(v):<5} median p={st.median(v):.3f}  min p={min(v):.2e}  p<2.5e-6: {sum((1 for x in v if x < 2.5e-06))}')
ap = [g for g in B if g.startswith('ENSG00000129354')]
if ap:
    p_ap, n_ap = B[ap[0]]
    peers = [p for n, p in pairs if abs(n - n_ap) <= 15]
    better = sum((1 for p in peers if p <= p_ap))
    print(f'\n=== AP1M2 vs same-size peers (n={n_ap}, +-15) ===')
    print(f'  peers: {len(peers)}   AP1M2 p={p_ap:.2e}')
    print(f'  peers with p <= AP1M2: {better}  ({better / len(peers) * 100:.2f}%)')
    print(f'  peer median p: {st.median(peers):.3f}')
gained = 0
gained_sig = 0
for g, (p, n) in B.items():
    if g in D:
        pd_, nd = D[g]
        if n > nd:
            gained += 1
            if p < 2.5e-06 and pd_ >= 2.5e-06:
                gained_sig += 1
print(f'\n=== rE2G gave extra variants to {gained} genes ===')
print(f'  of those, newly significant: {gained_sig}')
