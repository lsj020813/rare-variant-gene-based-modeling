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
f = [_required('COHORT_C_METADATA_FILE')][0]
z = zipfile.ZipFile(f)
x = z.read('xl/sharedStrings.xml').decode('utf-8', 'ignore')
ss = [html.unescape(re.sub('<[^>]+>', '', m)) for m in re.findall('<si>(.*?)</si>', x, re.S)]

def clean(s):
    return ' '.join(s.split())
conv = [i for i, s in enumerate(ss) if '전환' in s]
print(f'entries mentioning 전환: {len(conv)}')
for i in conv[:12]:
    print(f'  [{i}] {clean(ss[i])[:120]}')
print()
note = [i for i, s in enumerate(ss) if ('전환' in s and ('값' in s or '식' in s)) and len(s) > 40]
for i in note[:8]:
    print(f'  NOTE [{i}] {clean(ss[i])[:300]}')
print()
ori = [i for i, s in enumerate(ss) if 'Total cholesterol' in s or 'Total Cholesterol' in s]
for i in ori[:8]:
    print(f'  TCHL-label [{i}] {clean(ss[i])[:160]}')
