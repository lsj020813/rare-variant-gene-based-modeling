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
import collections, statistics as st, json
R = _configured('${PROJECT_ROOT}/work/ref/deductive')
STICKERS = ['splice_acceptor_variant', 'splice_donor_variant', 'stop_gained', 'frameshift_variant', 'stop_lost', 'start_lost', 'inframe_insertion', 'inframe_deletion', 'missense_variant', 'protein_altering_variant', 'splice_region_variant']
SSET = set(STICKERS)
pair = {}
cons = collections.Counter()
with open(f'{R}/vep_lof2.tsv') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 8 or f[7] != 'YES':
            continue
        hit = set(f[5].split(',')) & SSET
        if not hit:
            continue
        g = f[3].split('.')[0]
        k = f[0]
        v = pair.setdefault((g, k), [0] * 11)
        for cq in hit:
            v[STICKERS.index(cq)] = 1
            cons[cq] += 1
m = collections.Counter((g for g, _ in pair))
keys = {k for _, k in pair}
vals = sorted(m.values())
ge2 = {g for g, n in m.items() if n >= 2}
print(f'L0 variants (11-sticker, distinct keys): {len(keys):,}')
print(f'genes >=1: {len(m):,}   genes >=2: {len(ge2):,}')
print(f'm: median {st.median(vals)}  mean {st.mean(vals):.1f}  p90 {vals[int(len(vals) * 0.9)]}  max {vals[-1]}')
print('sticker counts:', {s: cons[s] for s in STICKERS})
with open(f'{R}/l0_variants.tsv', 'w') as out:
    out.write('gene\tkey37\t' + '\t'.join(STICKERS) + '\n')
    for (g, k), v in sorted(pair.items()):
        if g in ge2:
            out.write(f'{g}\t{k}\t' + '\t'.join(map(str, v)) + '\n')
n = sum((1 for _ in open(f'{R}/l0_variants.tsv'))) - 1
print(f'l0_variants.tsv rows (genes>=2 only): {n:,}')
print('L0SET_DONE')
