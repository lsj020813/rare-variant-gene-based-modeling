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
import glob, re, openpyxl
CB = _configured('${METADATA_DIR}')
path = [_required('COHORT_C_METADATA_FILE')][0]
wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
pat = re.compile(_configured('${COHORT_C_DISEASE_COLUMN_REGEX}'))
seen = {}
for ws in wb.worksheets:
    for row in ws.iter_rows(values_only=True):
        cells = [str(c).strip() for c in row if c is not None]
        for c in cells:
            if pat.match(c) and c not in seen:
                seen[c] = ' / '.join((x for x in cells if x != c))[:150]
wb.close()
for k in sorted(seen):
    print(k, '::', seen[k])
