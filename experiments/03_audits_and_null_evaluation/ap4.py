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
print('=== genes overlapping 10,680,000-10,810,000 ===')
with gzip.open(f'{R}/deductive/gencode.sorted.gtf.gz', 'rt') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.split('\t')
        if len(f) < 9 or f[2] != 'gene' or f[0] not in ('19', 'chr19'):
            continue
        s, e = (int(f[3]), int(f[4]))
        if e < 10680000 or s > 10810000:
            continue
        if 'gene_type "protein_coding"' not in f[8]:
            continue
        nm = f[8].split('gene_name "')[1].split('"')[0]
        gid = f[8].split('gene_id "')[1].split('"')[0]
        print(f'  {nm:<12} {gid:<22} {s:,}-{e:,}')

def load(fp):
    d = collections.defaultdict(set)
    for line in open(fp):
        f = line.split()
        if len(f) > 2 and f[1] == 'var':
            d[f[0]] = set(f[2:])
    return d
B = load(f'{R}/groupfiles_pilot/chr19.B_3kb_re2g.txt')
D = load(f'{R}/groupfiles_pilot/chr19.D_3kb_only.txt')
apB = next((v for k, v in B.items() if k.split('.')[0] == AP))
apD = next((v for k, v in D.items() if k.split('.')[0] == AP))
extra = apB - apD
owners = collections.Counter()
for g, vs in B.items():
    n = len(vs & extra)
    if n:
        owners[g] = n
print(f"\n=== who else owns AP1M2's {len(extra)} rE2G-only variants (arm B) ===")
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
for g, n in owners.most_common(10):
    print(f"  {name.get(g.split('.')[0], '?'):<12} {g:<22} shares {n}/{len(extra)}")
dnm = [g for g in B if name.get(g.split('.')[0]) == 'DNM2']
if dnm:
    dv = B[dnm[0]]
    print(f'\n=== AP1M2 vs DNM2 (arm B) ===')
    print(f'  AP1M2 {len(apB)}  DNM2 {len(dv)}  shared {len(apB & dv)}  jaccard {len(apB & dv) / len(apB | dv):.3f}')
    print(f"  of AP1M2's rE2G-only {len(extra)}: in DNM2 = {len(extra & dv)}")
