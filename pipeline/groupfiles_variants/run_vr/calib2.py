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
import subprocess, statistics as st
BCF = _configured('${BCFTOOLS}')
V = _configured('${GENOTYPE_DIR}/${GENOTYPE_FILENAME_PREFIX}22${GENOTYPE_FILENAME_SUFFIX}')
out = subprocess.run(f"{BCF} query -f '%POS\t%INFO/R2\t%INFO/AF\t[%DS ]\t[%GT ]\n' {V} 2>/dev/null | head -1500", shell=True, capture_output=True, text=True).stdout
rows = []
for line in out.strip().split('\n'):
    f = line.split('\t')
    if len(f) < 5:
        continue
    try:
        r2 = float(f[1])
        af = float(f[2])
    except ValueError:
        continue
    ds = [float(x) for x in f[3].split() if x not in ('.', '')]
    gt = f[4].split()
    ac = sum((g.count('1') for g in gt))
    rows.append((r2, af, sum(ds), ac, len(ds), gt[:3], ds[:3]))
print('=== is the >1000 band actually high-R2? ===')
for lo, hi, l in [(1, 10, 'DS 1-10'), (10, 20.5, 'DS 10-20.5'), (20.5, 100, 'DS 20-100'), (100, 1000, 'DS 100-1k'), (1000, 1000000000.0, 'DS >1k')]:
    s = [r for r in rows if lo <= r[2] < hi]
    if not s:
        continue
    print(f'  {l:<12} n={len(s):<5} median R2={st.median([r[0] for r in s]):.3f} median GT/DS={st.median([r[3] / r[2] for r in s if r[2] > 0.5]):.3f}')
print('\n=== sanity: a high-R2 common variant, raw values ===')
hi = [r for r in rows if r[0] > 0.95 and r[2] > 2000]
if hi:
    r = hi[0]
    print(f'  R2={r[0]:.3f} AF={r[1]:.4f} sum(DS)={r[2]:.1f} GTcount={r[3]} N={r[4]}')
    print(f'  expected AC from AF: {r[1] * 2 * r[4]:.1f}')
    print(f'  first GT: {r[5]}   first DS: {r[6]}')
else:
    print('  none with R2>0.95 and DS>2000 in this slice')
    top = sorted(rows, key=lambda x: -x[0])[:2]
    for r in top:
        print(f'  R2={r[0]:.3f} AF={r[1]:.4f} sum(DS)={r[2]:.1f} GT={r[3]} expAC={r[1] * 2 * r[4]:.1f} firstGT={r[5]} firstDS={r[6]}')
