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
import glob, gzip, collections
R = _configured('${PROJECT_ROOT}/work/ref')
AP = 'ENSG00000129354'
name = {}
with gzip.open(f'{R}/deductive/gencode.sorted.gtf.gz', 'rt') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.split('\t')
        if len(f) < 9 or f[2] != 'gene' or f[0] not in ('19', 'chr19'):
            continue
        if 'gene_name "' in f[8]:
            name[f[8].split('gene_id "')[1].split('"')[0].split('.')[0]] = f[8].split('gene_name "')[1].split('"')[0]

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
            sk = h.index('Pvalue_SKAT')
            bu = h.index('Pvalue_Burden')
            for line in fh:
                fl = line.rstrip('\n').split('\t')
                try:
                    out[fl[gi].split('.')[0]] = (float(fl[pi]), int(float(fl[ni])), float(fl[sk]), float(fl[bu]))
                except (ValueError, IndexError):
                    pass
    return out
B = read(sorted(glob.glob(f'{R}/saige_step2_pilot/tchl.chr19.B_3kb_re2g.B.part[0-9][0-9][0-9]')))
D = read([f'{R}/saige_step2_pilot/tchl.chr19.D_3kb_only'])
dsort = sorted(D.items(), key=lambda x: x[1][0])
rank = [i for i, (g, _) in enumerate(dsort, 1) if g == AP]
print(f'=== AP1M2 in arm D (3kb only, NO rE2G) ===')
print(f'  p={D[AP][0]:.2e}  n={D[AP][1]}  rank {rank[0]}/{len(D)}')
print(f'  arm D top 8:')
for g, (p, n, sk, bu) in dsort[:8]:
    mark = ' <-- AP1M2' if g == AP else ''
    print(f"    {name.get(g, '?'):<12} p={p:.2e} n={n}{mark}")
print(f'\n=== AP1M2 in arm B (with rE2G) ===')
print(f'  p={B[AP][0]:.2e}  n={B[AP][1]}   (threshold 2.5e-6)')
print(f'  SKAT {B[AP][2]:.2e} -> was {D[AP][2]:.2e}    Burden {B[AP][3]:.3f} -> was {D[AP][3]:.3f}')
print(f'  p improved {D[AP][0] / B[AP][0]:.1f}x by adding {B[AP][1] - D[AP][1]} variants')

def load(fp):
    d = {}
    for line in open(fp):
        f = line.split()
        if len(f) > 2 and f[1] == 'var':
            d[f[0].split('.')[0]] = set(f[2:])
    return d
gB = load(f'{R}/groupfiles_pilot/chr19.B_3kb_re2g.txt')
gD = load(f'{R}/groupfiles_pilot/chr19.D_3kb_only.txt')
extra = gB[AP] - gD[AP]
print(f"\n=== co-assigned genes (share AP1M2's {len(extra)} rE2G variants) ===")
rows = []
for g, vs in gB.items():
    s = len(vs & extra)
    if s >= 6 and g != AP and (g in B) and (g in D):
        rows.append((s, g, B[g][0], D[g][0]))
for s, g, pb, pd in sorted(rows, reverse=True):
    print(f"  {name.get(g, '?'):<12} shares {s:>2}/{len(extra)}   B p={pb:.2e}   D p={pd:.2e}")
