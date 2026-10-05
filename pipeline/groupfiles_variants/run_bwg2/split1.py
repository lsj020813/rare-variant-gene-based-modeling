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
import os
R = _configured('${PROJECT_ROOT}/work/ref')
SRC = f'{R}/groupfiles_bwg/chr1.B_3kb_re2g.txt'
OUT = f'{R}/groupfiles_bwg_chunks'
os.makedirs(OUT, exist_ok=True)
blocks = []
cur = []
genes = set()
slots = 0
for line in open(SRC):
    f = line.split()
    if f[1] == 'var':
        if cur:
            blocks.append(cur)
        cur = [line]
        genes.add(f[0])
        slots += len(f) - 2
    else:
        cur.append(line)
        assert f[2:] and set(f[2:]) == {'all'}, f'label not all: {f[:4]}'
if cur:
    blocks.append(cur)
assert len(blocks) == len(genes) == 2035, f'genes {len(blocks)}'
NCH = 18
per = (len(blocks) + NCH - 1) // NCH
w_genes = 0
w_slots = 0
for i in range(NCH):
    chunk = blocks[i * per:(i + 1) * per]
    if not chunk:
        break
    with open(f'{OUT}/chr1.part{i + 1:03d}.txt', 'w') as fh:
        for b in chunk:
            fh.write(''.join(b))
            w_genes += 1
            w_slots += len(b[0].split()) - 2
assert w_genes == 2035, f'union genes {w_genes}'
assert w_slots == slots, f'union slots {w_slots} vs {slots}'
n = len([f for f in os.listdir(OUT) if f.startswith('chr1.part')])
print(f'chunks {n}  genes {w_genes}/2035  slots {w_slots}/{slots}  label=all verified')
print('CHR1_SPLIT_OK')
