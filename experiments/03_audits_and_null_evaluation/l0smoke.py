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
import subprocess, numpy as np, collections
R = _configured('${PROJECT_ROOT}/work/ref')
BCF = _configured('${BCFTOOLS}')
rows = [l.rstrip('\n').split('\t') for l in open(f'{R}/deductive/l0_variants.tsv')][1:]
ch21 = [(g, k) for g, k, *st in rows if k.split(':')[0] == 'chr21']
keys21 = {k for _, k in ch21}
print(f'chr21 L0 variants: {len(keys21):,} ({len(ch21):,} gene-variant rows)')
pos = sorted({int(k.split(':')[1]) for k in keys21})
regions = ','.join((f'21:{p}-{p}' for p in pos[:50]))
out = subprocess.run(f"{BCF} query -r {regions} -f '%CHROM:%POS:%REF:%ALT[\t%DS]\n' {R}/band_vcf/chr21.band.vcf.gz", shell=True, capture_output=True, text=True)
got = 0
nsamp = None
matched = 0
for line in out.stdout.splitlines():
    f = line.split('\t')
    got += 1
    if nsamp is None:
        nsamp = len(f) - 1
    if 'chr21:' + f[0].split(':', 1)[1] in keys21 or f[0] in keys21:
        matched += 1
print(f'fetched rows: {got}  samples/row: {nsamp}  key-matched: {matched}')
print('first key from vcf:', out.stdout.split('\t')[0] if out.stdout else 'EMPTY')
print('first key from l0 table:', sorted(keys21)[0])
assert got > 0 and nsamp == N_SAMPLES, f'GATE FAIL: rows={got} nsamp={nsamp}'
print('SMOKE_EXTRACT_OK')
