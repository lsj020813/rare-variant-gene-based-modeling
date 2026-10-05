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
import csv, os, re
B = _configured('${PHENO_DIR}')
src = open(_configured('${PROJECT_ROOT}/work/run_pheno2/build.py'), encoding='utf-8', errors='ignore').read()
for m in re.finditer('(HT|DM|LP|DIHT|PDHT|TRAIT|DIS)\\w*', src):
    pass
import itertools
lines = [l for l in src.splitlines() if re.search('HT|DM|LP|DISEASE|glob', l)]
print('=== v2 build lines mentioning traits/disease ===')
for l in lines[:22]:
    print('  ', l.strip()[:150])
print()
print('=== actual disease-table columns per cohort ===')
for rel in [_configured('${COHORT_A_DISEASE_FILE}'), _configured('${COHORT_C_DISEASE_FILE}'), _configured('${COHORT_B_DISEASE_FILE}')]:
    p = _source_path(B, rel)
    if not os.path.exists(p):
        print(f'  {rel}: MISSING')
        continue
    with open(p, encoding='utf-8', errors='ignore') as fh:
        h = fh.readline().rstrip('\n').split('\t')
    hit = [c for c in h if re.search(_required('DISEASE_COLUMN_REGEX'), c, re.I)]
    print(f'  [{rel}] {len(h)} cols -> {hit[:12]}')
