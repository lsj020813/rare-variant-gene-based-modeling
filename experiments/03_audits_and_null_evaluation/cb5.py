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
import glob, zipfile, re, html
D = _configured('${METADATA_DIR}')
for metadata_path, label in [(_required('COHORT_A_METADATA_FILE'), 'CT'), (_required('COHORT_B_METADATA_FILE'), 'NC')]:
    f = [metadata_path]
    if not f:
        print(f'[{label}] metadata NOT FOUND')
        continue
    z = zipfile.ZipFile(f[0])
    ss = [html.unescape(re.sub('<[^>]+>', '', m)) for m in re.findall('<si>(.*?)</si>', z.read('xl/sharedStrings.xml').decode('utf-8', 'ignore'), re.S)]
    print(f"\n=== {label}: {f[0].split('/')[-1]} ({len(ss)} strings) ===")
    conv = [s for s in ss if '전환' in s]
    print(f'  entries mentioning 전환(conversion): {len(conv)}')
    for s in conv[:3]:
        print('   ', ' '.join(s.split())[:150])
    tc = [s for s in ss if 'TCHL' in s.upper() or 'Total cholesterol' in s or 'Total Cholesterol' in s]
    print(f'  TCHL entries: {len(tc)}')
    for s in tc[:4]:
        print('   ', ' '.join(s.split())[:200])
    dev = [s for s in ss if '장비' in s]
    print(f'  entries mentioning 장비(instrument): {len(dev)}')
    for s in dev[:2]:
        print('   ', ' '.join(s.split())[:180])
