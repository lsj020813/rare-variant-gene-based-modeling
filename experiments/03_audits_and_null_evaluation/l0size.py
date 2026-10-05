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
import glob, collections
R = _configured('${PROJECT_ROOT}/work/ref')
genes = collections.Counter()
keys = set()
cons = collections.Counter()
for fp in glob.glob(f'{R}/deductive/chr*.a.tsv'):
    for line in open(fp):
        f = line.rstrip('\n').split('\t')
        if len(f) < 4:
            continue
        k = f[0]
        g = f[2].split('.')[0]
        cq = f[3]
        keys.add(k)
        genes[g] += 0
        cons[cq] += 1
pair = set()
for fp in glob.glob(f'{R}/deductive/chr*.a.tsv'):
    for line in open(fp):
        f = line.rstrip('\n').split('\t')
        if len(f) < 4:
            continue
        pair.add((f[2].split('.')[0], f[0]))
m = collections.Counter((g for g, _ in pair))
import statistics as st
vals = sorted(m.values())
print(f'coding variants (distinct keys): {len(keys):,}')
print(f'genes with >=1 coding band variant: {len(m):,}')
print(f'genes with >=2: {sum((1 for v in vals if v >= 2)):,}')
print(f'm per gene: median {st.median(vals)}  mean {st.mean(vals):.1f}  p90 {vals[int(len(vals) * 0.9)]}  max {vals[-1]}')
print('top consequences:', cons.most_common(8))
