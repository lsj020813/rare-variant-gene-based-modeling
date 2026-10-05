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
import glob, zipfile, re, html, csv, collections
D = _configured('${METADATA_DIR}')
f = [_required('COHORT_A_METADATA_FILE')][0]
z = zipfile.ZipFile(f)
ss = [html.unescape(re.sub('<[^>]+>', '', m)) for m in re.findall('<si>(.*?)</si>', z.read('xl/sharedStrings.xml').decode('utf-8', 'ignore'), re.S)]
i = [k for k, s in enumerate(ss) if s.strip() == _required('COHORT_A_TCHL_INSTRUMENT_COLUMN')]
print('instrument column metadata context:')
for k in i:
    for j in range(max(0, k - 2), min(len(ss), k + 4)):
        print(f"  [{j}] {' '.join(ss[j].split())[:220]}")
print()
B = _configured('${PHENO_DIR}')
with open(_source_path(B, _configured('${COHORT_A_BIOCHEM_FILE}')), encoding='utf-8', errors='ignore') as fh:
    rd = csv.reader(fh, delimiter='\t')
    hdr = next(rd)
    ii = hdr.index(_required('COHORT_A_TCHL_INSTRUMENT_COLUMN'))
    tt = hdr.index(_required('COHORT_A_TCHL_COLUMN'))
    cnt = collections.Counter()
    byinst = collections.defaultdict(list)
    for r in rd:
        if len(r) <= max(ii, tt):
            continue
        v = r[ii].strip()
        cnt[v] += 1
        try:
            byinst[v].append(float(r[tt]))
        except ValueError:
            pass
import statistics as st
print('instrument column value counts and per-value TCHL mean:')
for v, n in cnt.most_common(12):
    m = st.mean(byinst[v]) if byinst[v] else None
    print(f"  '{v}': n={n:,}  mean TCHL={m:.1f}" if m else f"  '{v}': n={n:,}")
