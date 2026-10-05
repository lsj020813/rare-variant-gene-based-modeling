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
import glob, gzip, re, collections
R = _configured('${PROJECT_ROOT}/work/ref')
THR = 2.5e-06
sym = {}
with gzip.open(f'{R}/deductive/gencode.sorted.gtf.gz', 'rt') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.split('\t', 9)
        if len(f) < 9 or f[2] != 'gene':
            continue
        a = f[8]
        if 'gene_name "' in a:
            sym[a.split('gene_id "')[1].split('"')[0].split('.')[0]] = a.split('gene_name "')[1].split('"')[0]
CHRS = ['15', '16', '17', '18', '19', '20', '21', '22']

def read_arm(paths):
    out = {}
    for fp in paths:
        with open(fp) as fh:
            h = fh.readline().rstrip('\n').split('\t')
            if 'Pvalue' not in h:
                continue
            gi = h.index('Region')
            pi = h.index('Pvalue')
            for line in fh:
                f = line.rstrip('\n').split('\t')
                try:
                    out[f[gi].split('.')[0]] = float(f[pi])
                except (ValueError, IndexError):
                    pass
    return out
part_re = re.compile('\\.part\\d{3}$')
print(f"{'trait':<6}{'chr':<6}{'B sig':>6}{'v4 sig':>7}   B-arm significant genes")
tot_b = collections.Counter()
tot_v = collections.Counter()
for T in ('tchl', 'htn', 'dm', 'lip'):
    for CH in CHRS:
        b = read_arm([p for p in [f'{R}/saige_step2_bwg/{T}.chr{CH}'] if glob.glob(p)])
        v = read_arm([p for p in glob.glob(f'{R}/saige_step2_v4/{T}.chr{CH}.part*') if part_re.search(p)])
        bs = [g for g, p in b.items() if p < THR]
        vs = [g for g, p in v.items() if p < THR]
        tot_b[T] += len(bs)
        tot_v[T] += len(vs)
        if bs or vs:
            names = ','.join((f'{sym.get(g, g)}({b[g]:.0e})' for g in sorted(bs, key=lambda x: b[x])))
            print(f'{T:<6}chr{CH:<4}{len(bs):>5}{len(vs):>7}   {names[:100]}')
print()
print('totals (chr15-22):')
for T in ('tchl', 'htn', 'dm', 'lip'):
    print(f'  {T}: B-arm {tot_b[T]}   v4(100kb) {tot_v[T]}')
