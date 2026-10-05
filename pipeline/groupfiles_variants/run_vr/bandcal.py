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
B = _configured('${PROJECT_ROOT}/work/ref/band_vcf/chr22.band.vcf.gz')
hdr = subprocess.run(f"{BCF} view -h {B} | grep '^##FORMAT'", shell=True, capture_output=True, text=True).stdout
print('=== band FORMAT fields ===')
print(hdr.strip()[:300])
out = subprocess.run(f"{BCF} query -f '%INFO/R2\t%INFO/MAF\t[%DS ]\t[%GT ]\n' {B} 2>/dev/null | head -4000", shell=True, capture_output=True, text=True).stdout
rows = []
for line in out.strip().split('\n'):
    f = line.split('\t')
    if len(f) < 4:
        continue
    try:
        r2 = float(f[0])
        maf = float(f[1])
    except ValueError:
        continue
    ds = [float(x) for x in f[2].split() if x not in ('.', '')]
    gt = f[3].split()
    rows.append((r2, maf, sum(ds), sum((g.count('1') for g in gt))))
print(f'\nband variants sampled: {len(rows)}')
if rows:
    print(f'R2: min {min((r[0] for r in rows)):.3f}  median {st.median([r[0] for r in rows]):.3f}  max {max((r[0] for r in rows)):.3f}')
    print('\n=== GT/DS agreement by R2 ===')
    for lo, hi in [(0.3, 0.5), (0.5, 0.7), (0.7, 0.9), (0.9, 1.01)]:
        s = [r for r in rows if lo <= r[0] < hi and r[2] > 0.5]
        if not s:
            continue
        rt = [r[3] / r[2] for r in s]
        print(f'  R2 {lo:.1f}-{hi:.1f}: n={len(s):<5} median GT/DS={st.median(rt):.3f}')
    n_lo = sum((1 for r in rows if r[0] < 0.5))
    print(f'\nband share with R2<0.5: {n_lo / len(rows) * 100:.1f}%')
    allr = [r[3] / r[2] for r in rows if r[2] > 0.5]
    print(f'band overall median GT/DS: {st.median(allr):.3f}')
