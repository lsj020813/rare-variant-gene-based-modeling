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
vcf = {}
out = subprocess.run(f"{BCF} query -f '%POS\t%INFO/MAF\t%INFO/R2\n' {W}/geno_smoke/chr22.vcf.gz", shell=True, capture_output=True, text=True).stdout
for l in out.strip().split('\n'):
    f = l.split('\t')
    if len(f) == 3:
        vcf[f[0]] = (float(f[1]), float(f[2]))
bim = [l.split('\t')[3] for l in open(f'{W}/geno_smoke/chr22.bim')]
ac = []
with open(f'{W}/geno_smoke/chr22.f.acount') as fh:
    fh.readline()
    for l in fh:
        f = l.split('\t')
        ac.append((float(f[5]), float(f[6])))
print(f'bim {len(bim)}  acount {len(ac)}')
print(f"{'POS':>10} {'INFO/MAF':>10} {'MAF*2N':>8} {'hardcall_MAC':>13} {'R2':>6}")
mism = 0
for i, pos in enumerate(bim[:10]):
    if pos in vcf and i < len(ac):
        maf, r2 = vcf[pos]
        alt, obs = ac[i]
        print(f'{pos:>10} {maf:>10.3e} {maf * 174860:>8.1f} {alt:>13.0f} {r2:>6.2f}')
inb = sum((1 for a, o in ac if 10 < a <= 20.5))
print(f'\nhard-call MAC 10-20.5: {inb}/{len(ac)}')
print(f'hard-call MAC <=10   : {sum((1 for a, o in ac if a <= 10))}/{len(ac)}')
