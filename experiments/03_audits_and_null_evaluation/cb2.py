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
print('file:', f.split('/')[-1])
z = zipfile.ZipFile(f)
ss = []
if 'xl/sharedStrings.xml' in z.namelist():
    x = z.read('xl/sharedStrings.xml').decode('utf-8', 'ignore')
    ss = [html.unescape(re.sub('<[^>]+>', '', m)) for m in re.findall('<si>(.*?)</si>', x, re.S)]
print('shared strings:', len(ss))
idx = [i for i, s in enumerate(ss) if 'TCHL' in s.upper()]
print('entries containing TCHL:', len(idx))
for i in idx[:20]:
    ctx = [ss[j] for j in range(max(0, i - 3), min(len(ss), i + 5))]
    print(f'  [{i}] ' + ' | '.join((c.replace(chr(10), ' ')[:70] for c in ctx)))
