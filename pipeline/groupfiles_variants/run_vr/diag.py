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
import os as _os
N_SAMPLES = _required_int('N_SAMPLES')
import subprocess, csv
W = _configured('${PROJECT_ROOT}/work/run_vr')
P2 = _configured('${PLINK2}')
sites = {}
for l in open(f'{W}/chunks/chr22.rare.tsv'):
    f = l.split('\t')
    sites[f[0], f[1]] = None
BCF = _configured('${BCFTOOLS}')
out = subprocess.run(f"{BCF} query -f '%CHROM\t%POS\t%INFO/MAF\t%INFO/AF\t%INFO/R2\n' {W}/geno_smoke/chr22.vcf.gz", shell=True, capture_output=True, text=True).stdout
vcf = {}
for l in out.strip().split('\n'):
    f = l.split('\t')
    if len(f) >= 3:
        vcf[f[0], f[1]] = (float(f[2]), float(f[3]) if f[3] not in ('', '.') else None)
print(f'VCF sites: {len(vcf)}')
ac = {}
with open(f'{W}/geno_smoke/chr22.f.acount') as fh:
    h = fh.readline().split()
    ci = h.index('ALT_CTS') if 'ALT_CTS' in h else 5
    oi = h.index('OBS_CT') if 'OBS_CT' in h else 6
    pi = h.index('POS') if 'POS' in h else 1
    for l in fh:
        f = l.split()
        ac[f[pi]] = (float(f[ci]), float(f[oi]))
print(f'acount rows: {len(ac)}')
print(f"\n{'POS':>10} {'INFO/MAF':>10} {'expect MAC':>11} {'plink ALT_CTS':>14} {'OBS_CT':>8}")
n = 0
for (ch, pos), (maf, af) in list(vcf.items()):
    if pos in ac and n < 8:
        alt, obs = ac[pos]
        print(f'{pos:>10} {maf:>10.3e} {maf * 2 * N_SAMPLES:>11.1f} {alt:>14.0f} {obs:>8.0f}')
        n += 1
