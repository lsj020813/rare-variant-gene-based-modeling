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
import collections, statistics as st
R = _configured('${PROJECT_ROOT}/work/ref/deductive')
CODING = {'missense_variant', 'synonymous_variant', 'stop_gained', 'stop_lost', 'start_lost', 'frameshift_variant', 'inframe_insertion', 'inframe_deletion', 'splice_donor_variant', 'splice_acceptor_variant', 'stop_retained_variant', 'protein_altering_variant', 'incomplete_terminal_codon_variant', 'coding_sequence_variant'}
pair = set()
cons = collections.Counter()
with open(f'{R}/vep_lof2.tsv') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 8 or f[7] != 'YES':
            continue
        cqs = set(f[5].split(','))
        hit = cqs & CODING
        if not hit:
            continue
        g = f[3].split('.')[0]
        k = f[0]
        pair.add((g, k))
        for cq in hit:
            cons[cq] += 1
m = collections.Counter((g for g, _ in pair))
keys = {k for _, k in pair}
vals = sorted(m.values())
print(f'L0 coding variants (distinct): {len(keys):,}')
print(f'genes >=1: {len(m):,}   genes >=2: {sum((1 for v in vals if v >= 2)):,}')
if vals:
    print(f'm: median {st.median(vals)}  mean {st.mean(vals):.1f}  p90 {vals[int(len(vals) * 0.9)]}  max {vals[-1]}')
print('consequences:', cons.most_common(10))
