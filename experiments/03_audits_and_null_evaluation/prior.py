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
import csv, collections, glob, os
D = _configured('${SAIGE_PHENOTYPE_DIR}')
for f in sorted(glob.glob(f'{D}/*.saige.tsv')):
    with open(f, encoding='utf-8', errors='ignore') as fh:
        rd = csv.reader(fh, delimiter='\t')
        hdr = next(rd)
        rows = [r for r in rd if r]
    print(f'\n=== {os.path.basename(f)} ===')
    print('  columns:', hdr)
    print('  n:', f'{len(rows):,}')
    for col in ('CT', 'NC', 'AS', 'cohort', 'batch', 'INST', 'inst'):
        if col in hdr:
            i = hdr.index(col)
            print(f'  {col} value counts:', dict(collections.Counter((r[i] for r in rows)).most_common(6)))
