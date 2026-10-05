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
print('=== protein-coding genes 10.50-11.30 Mb (chr19) ===')
loc = {}
with gzip.open(f'{R}/deductive/gencode.sorted.gtf.gz', 'rt') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.split('\t')
        if len(f) < 9 or f[2] != 'gene' or f[0] not in ('19', 'chr19'):
            continue
        if 'gene_type "protein_coding"' not in f[8]:
            continue
        s, e = (int(f[3]), int(f[4]))
        if e < 10500000 or s > 11300000:
            continue
        nm = f[8].split('gene_name "')[1].split('"')[0]
        gid = f[8].split('gene_id "')[1].split('"')[0]
        loc[gid.split('.')[0]] = (nm, s, e)
for g, (nm, s, e) in sorted(loc.items(), key=lambda x: x[1][1]):
    star = ' ***' if nm in ('LDLR', 'SMARCA4', 'CARM1', 'AP1M2', 'SLC44A2', 'DNM2') else ''
    print(f'  {nm:<12} {s / 1000000.0:>7.3f}-{e / 1000000.0:.3f} Mb{star}')

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
                    out[fl[gi].split('.')[0]] = (float(fl[pi]), int(float(fl[ni])))
                except (ValueError, IndexError):
                    pass
    return out
A = read(sorted(glob.glob(f'{R}/saige_step2_v4/tchl.chr19.part[0-9][0-9][0-9]')))
B = read(sorted(glob.glob(f'{R}/saige_step2_pilot/tchl.chr19.B_3kb_re2g.B.part[0-9][0-9][0-9]')))
C = read([f'{R}/saige_step2_pilot/tchl.chr19.C_nearest'])
D = read([f'{R}/saige_step2_pilot/tchl.chr19.D_3kb_only'])
print('\n=== p-values in this window (blank = not tested) ===')
print(f"  {'gene':<12}{'A 100kb':>12}{'B 3kb+rE2G':>13}{'C nearest':>12}{'D 3kb':>12}")
for g, (nm, s, e) in sorted(loc.items(), key=lambda x: x[1][1]):

    def f(d):
        return f'{d[g][0]:.1e}' if g in d else '-'
    row = [f(A), f(B), f(C), f(D)]
    if all((r == '-' for r in row)):
        continue
    best = min([d[g][0] for d in (A, B, C, D) if g in d])
    mark = ' <<<' if best < 2.5e-06 else ''
    print(f'  {nm:<12}{row[0]:>12}{row[1]:>13}{row[2]:>12}{row[3]:>12}{mark}')
ldlr = [g for g, (nm, _, _) in loc.items() if nm == 'LDLR']
if ldlr:
    L = ldlr[0]
    print(f'\n=== LDLR ({L}) ===')
    for tag, d in (('A', A), ('B', B), ('C', C), ('D', D)):
        print(f'  {tag}: p={d[L][0]:.2e} n={d[L][1]}' if L in d else f'  {tag}: NOT TESTED')

    def load(fp):
        out = {}
        for line in open(fp):
            f = line.split()
            if len(f) > 2 and f[1] == 'var':
                out[f[0].split('.')[0]] = set(f[2:])
        return out
    gB = load(f'{R}/groupfiles_pilot/chr19.B_3kb_re2g.txt')
    if L in gB and 'ENSG00000129354' in gB:
        a = gB['ENSG00000129354']
        l = gB[L]
        print(f'  AP1M2 {len(a)} vars, LDLR {len(l)} vars, shared {len(a & l)}')
