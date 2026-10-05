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
cmd = f"{BCF} query -f '%POS\t%INFO/R2\t%INFO/IMPUTED\t[%DS ]\t[%GT ]\n' {V} 2>/dev/null | head -3000"
out = subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout
rows = []
for line in out.strip().split('\n'):
    f = line.split('\t')
    if len(f) < 5:
        continue
    try:
        r2 = float(f[1])
    except ValueError:
        r2 = None
    imp = f[2] != '.'
    ds = [float(x) for x in f[3].split() if x not in ('.', '')]
    gt = f[4].split()
    ac = sum((g.count('1') for g in gt))
    rows.append((f[0], r2, imp, sum(ds), ac, len(ds)))
print(f'variants parsed: {len(rows)}')
gtd = [r for r in rows if not r[2]]
imp = [r for r in rows if r[2]]
print(f'  genotyped (truth): {len(gtd)}   imputed: {len(imp)}')

def report(tag, sub, lo=None, hi=None):
    s = [r for r in sub if lo is None or (r[3] is not None and lo <= r[3] < hi)]
    if not s:
        return
    ratios = [r[4] / r[3] for r in s if r[3] > 0.5]
    if not ratios:
        return
    print(f'  {tag:<26} n={len(s):<5} median GT/DS ratio = {st.median(ratios):.3f}')
print('\n=== GENOTYPED variants: does GT count equal DS sum? ===')
report('all genotyped', gtd)
print('\n=== IMPUTED variants, by dosage-MAC band ===')
for lo, hi, lbl in [(1, 10, 'DS sum 1-10'), (10, 20.5, 'DS sum 10-20.5'), (20.5, 100, 'DS sum 20-100'), (100, 1000, 'DS sum 100-1000'), (1000, 1000000000.0, 'DS sum >1000')]:
    report(lbl, imp, lo, hi)
print('\n=== IMPUTED, by R2 (DS sum 10-200 only) ===')
band = [r for r in imp if r[3] and 10 <= r[3] < 200 and (r[1] is not None)]
for lo, hi in [(0.0, 0.3), (0.3, 0.5), (0.5, 0.7), (0.7, 0.9), (0.9, 1.01)]:
    s = [r for r in band if lo <= r[1] < hi]
    if not s:
        continue
    ratios = [r[4] / r[3] for r in s if r[3] > 0.5]
    print(f'  R2 {lo:.1f}-{hi:.1f}: n={len(s):<5} median GT/DS = {st.median(ratios):.3f}')
