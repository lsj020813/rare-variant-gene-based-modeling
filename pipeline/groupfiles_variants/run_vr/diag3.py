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
import subprocess
W = _configured('${PROJECT_ROOT}/work/run_vr')
BCF = _configured('${BCFTOOLS}')
ds = subprocess.run(f"{BCF} query -f '[%DS ]\n' {W}/geno_smoke/chr22.vcf.gz", shell=True, capture_output=True, text=True).stdout.split('\n')[0].split()
gt = subprocess.run(f"{BCF} query -f '[%GT ]\n' {W}/geno_smoke/chr22.vcf.gz", shell=True, capture_output=True, text=True).stdout.split('\n')[0].split()
v = [float(x) for x in ds if x not in ('.', '')]
print(f'samples: {len(v)}   sum(DS) = {sum(v):.2f}  -> dosage-based MAC')
b = [0] * 5
for x in v:
    if x < 0.05:
        b[0] += 1
    elif x < 0.5:
        b[1] += 1
    elif x < 0.95:
        b[2] += 1
    elif x < 1.5:
        b[3] += 1
    else:
        b[4] += 1
print(f'  DS<0.05: {b[0]:,}   0.05-0.5: {b[1]}   0.5-0.95: {b[2]}   0.95-1.5: {b[3]}   >=1.5: {b[4]}')
alt = sum((g.count('1') for g in gt))
print(f'  GT alt alleles (hard call in the FILE): {alt}')
keep = sum((1 for x in v if abs(x - round(x)) <= 0.2))
print(f'  within plink default hard-call window (|x-round(x)|<=0.2): {keep:,}/{len(v):,}')
print(f'  would be MISSING at default: {len(v) - keep:,}')
