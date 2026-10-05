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
import glob, subprocess
BCF = _configured('${BCFTOOLS}')
G = _configured('${PROJECT_ROOT}/work/ref/groupfiles_chunks')
P = _configured('${PROJECT_ROOT}/work/ref/groupfiles_pilot')

def pool(pat):
    s = set()
    for fp in glob.glob(pat):
        for line in open(fp):
            f = line.split()
            if len(f) > 2 and f[1] == 'var':
                s.update(f[2:])
    return s
A = pool(f'{G}/chr19.part*.txt')
print(f'A (baseline, running) distinct variants: {len(A):,}')
arms = {}
for arm in ('B_3kb_re2g', 'C_nearest', 'D_3kb_only'):
    s = pool(f'{P}/chr19.{arm}.txt')
    arms[arm] = s
    out = len(s - A)
    print(f'{arm:<14} distinct={len(s):<8} outside_A={out:<6} subset_of_A={out == 0}')
ids = subprocess.run(_configured(f"{BCF} query -f '%ID\\n' ${{PROJECT_ROOT}}/work/ref/band_vcf/chr19.band.vcf.gz"), shell=True, capture_output=True, text=True).stdout.split()
band = set((i.strip() for i in ids if i.strip() and i != '.'))
print(f'\nband VCF chr19 IDs: {len(band):,}')
print(f'  A outside band: {len(A - band):,}')
for arm, s in arms.items():
    print(f'  {arm:<14} outside band: {len(s - band):,}')
print(f'\nsample of band ID: {next(iter(band))}')
print(f'sample of A key   : {next(iter(A))}')
