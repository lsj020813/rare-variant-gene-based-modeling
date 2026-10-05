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
import gzip, subprocess
R = _configured('${PROJECT_ROOT}/work/ref')
gf = f'{R}/groupfiles_bwg/chr19.B_3kb_re2g.txt'
vars_g = []
for line in open(gf):
    f = line.split()
    if f[0].split('.')[0] == 'ENSG00000129354' and f[1] == 'var':
        vars_g = f[2:]
        break
pos = [int(v.split(':')[1]) for v in vars_g]
print(f'AP1M2 variants n={len(pos)}  key-pos range {min(pos):,} - {max(pos):,}')
print(f'Configured b38 gene span: {_required_int("GENE_START_B38"):,} - {_required_int("TSS_POSITION"):,}')
BCF = _configured('${BCFTOOLS}')
k = vars_g[0]
out = subprocess.run(f"{BCF} query -f '%CHROM\t%POS\t%ID\n' {R}/lift38_keyed/chr19.keyed38.vcf.gz | grep -m3 -F '{k}'", shell=True, capture_output=True, text=True)
print('keyed map lookup for first key:', out.stdout.strip() or '(not found)')
