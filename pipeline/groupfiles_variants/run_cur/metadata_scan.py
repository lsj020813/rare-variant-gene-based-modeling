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
import glob, json, re
import openpyxl
CB = _configured('${METADATA_DIR}')
pat = re.compile(_required('METADATA_COLUMN_REGEX'))
out = {}
for path in [_required('COHORT_A_METADATA_FILE'), _required('COHORT_B_METADATA_FILE'), _required('COHORT_C_METADATA_FILE')]:
    name = path.split('/')[-1]
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    for ws in wb.worksheets:
        for row in ws.iter_rows(values_only=True):
            cells = [str(c).strip() for c in row if c is not None]
            if not cells:
                continue
            for c in cells:
                if pat.match(c):
                    key = f"{name.split('_')[0]}|{c}"
                    txt = ' / '.join((x for x in cells if x != c))[:400]
                    out.setdefault(key, txt)
    wb.close()
json.dump(out, open(_configured('${PROJECT_ROOT}/work/run_cur/metadata_hits.json'), 'w'), ensure_ascii=False, indent=1)
print(f'hits: {len(out)}')
for k in sorted(out):
    print(f'{k}  ::  {out[k][:200]}')
