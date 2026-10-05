#!/usr/bin/env python3
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
import glob, openpyxl
CB = _configured('${METADATA_DIR}')
p = [_required('COHORT_C_METADATA_FILE')][0]
wb = openpyxl.load_workbook(p, read_only=True, data_only=True)
ws = wb[_required('METADATA_SHEET')]
hits = 0
for i, row in enumerate(ws.iter_rows(values_only=True)):
    cells = [str(c).strip() for c in row if c is not None]
    joined = ' | '.join(cells)
    if _required('METADATA_SEARCH_TERM') in joined:
        print(f'[r{i}]', joined[:220])
        hits += 1
        if hits >= 8:
            break
wb.close()
print('done, hits:', hits)
