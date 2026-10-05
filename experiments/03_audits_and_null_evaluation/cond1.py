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
import subprocess, os, glob
BCF = _configured('${BCFTOOLS}')
BAND = _configured('${PROJECT_ROOT}/work/ref/band_vcf')
RAW = _configured('${GENOTYPE_DIR}')
REG = {'19': [(45411941, 45412079, 'APOE e4 rs429358 area'), (11202306, 11202400, 'LDLR rs6511720 area')], '1': [(109817590, 109818530, 'CELSR2/SORT1 rs646776 area')], '16': [(56993324, 56993400, 'CETP rs3764261 area')], '5': [(74650000, 74660000, 'HMGCR-ish / chr5 locus')]}
print('=== raw (GRCh37) genotype files present? ===')
for ch in ('1', '5', '16', '19'):
    hits = glob.glob(_configured(f'{RAW}/${{GENOTYPE_FILENAME_PREFIX}}{ch}.*.vcf.gz'))
    print(f"  chr{ch}: {len(hits)} file(s)  {(os.path.basename(hits[0]) if hits else '-')}")
print('\n=== are common variants (MAF>5%) present at these loci in RAW? ===')
for ch, regs in REG.items():
    f = glob.glob(_configured(f'{RAW}/${{GENOTYPE_FILENAME_PREFIX}}{ch}.*.vcf.gz'))
    if not f:
        print(f'  chr{ch}: no file')
        continue
    for s, e, lab in regs:
        out = subprocess.run(f"{BCF} query -r {ch}:{s}-{e} -f '%POS\\t%REF\\t%ALT\\t%INFO/MAF\\t%INFO/R2\\n' {f[0]} 2>/dev/null | head -4", shell=True, capture_output=True, text=True).stdout.strip()
        n = subprocess.run(f"{BCF} query -r {ch}:{s}-{e} -f '%POS\\n' {f[0]} 2>/dev/null | wc -l", shell=True, capture_output=True, text=True).stdout.strip()
        print(f'  chr{ch} {lab}: {n} variants')
        for l in out.split('\n')[:3]:
            if l:
                print(f'      {l}')
print('\n=== same window in our BAND file (should be near-empty: common variants excluded) ===')
for ch, regs in REG.items():
    bf = f'{BAND}/chr{ch}.band.vcf.gz'
    if not os.path.exists(bf):
        continue
    for s, e, lab in regs[:1]:
        n = subprocess.run(f"{BCF} query -r {ch}:{s}-{e} -f '%POS\\n' {bf} 2>/dev/null | wc -l", shell=True, capture_output=True, text=True).stdout.strip()
        print(f'  chr{ch} {lab}: {n} band variants')
