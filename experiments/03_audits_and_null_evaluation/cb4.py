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
import glob, zipfile, re, html, csv, statistics as st
D = _configured('${METADATA_DIR}')
f = [_required('COHORT_C_METADATA_FILE')][0]
z = zipfile.ZipFile(f)
ss = [html.unescape(re.sub('<[^>]+>', '', m)) for m in re.findall('<si>(.*?)</si>', z.read('xl/sharedStrings.xml').decode('utf-8', 'ignore'), re.S)]
print('=== TCHL 전환값 full text ===')
print(' '.join(ss[int(_required('METADATA_TRANSFORM_LABEL_INDEX'))].split()))
print()
print('=== TCHL original label full text ===')
print(' '.join(ss[int(_required('METADATA_ORIGINAL_LABEL_INDEX'))].split()))
print()
B = _configured('${PHENO_DIR}')
p = _source_path(B, _configured('${COHORT_C_BIOCHEM_FILE}'))
with open(p, encoding='utf-8', errors='ignore') as fh:
    rd = csv.reader(fh, delimiter='\t')
    hdr = next(rd)
    io, it = (hdr.index(_required('COHORT_C_TCHL_ORIGINAL_COLUMN')), hdr.index(_required('COHORT_C_TCHL_COLUMN')))
    pairs = []
    for r in rd:
        if len(r) <= max(io, it):
            continue
        try:
            o, t = (float(r[io]), float(r[it]))
        except ValueError:
            continue
        pairs.append((o, t))
print(f'paired rows: {len(pairs):,}')
same = sum((1 for o, t in pairs if abs(o - t) < 1e-09))
print(f'  ORI == TR (unchanged): {same:,} ({same / len(pairs) * 100:.1f}%)')
print(f'  ORI != TR (converted): {len(pairs) - same:,} ({(len(pairs) - same) / len(pairs) * 100:.1f}%)')
diff = [(o, t) for o, t in pairs if abs(o - t) >= 1e-09]
if diff:
    n = len(diff)
    mo = st.mean((o for o, _ in diff))
    mt = st.mean((t for _, t in diff))
    b = sum(((o - mo) * (t - mt) for o, t in diff)) / sum(((o - mo) ** 2 for o, _ in diff))
    a = mt - b * mo
    resid = [abs(t - (a + b * o)) for o, t in diff]
    print(f'  fitted on converted subset: TR = {a:.4f} + {b:.4f} * ORI')
    print(f'  max |residual|: {max(resid):.6f}  (exact linear conversion if ~0)')
    print(f'  ORI mean {mo:.2f} -> TR mean {mt:.2f}')
