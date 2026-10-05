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
import subprocess, os, collections
R = _configured('${PROJECT_ROOT}/work/ref')
BCF = _configured('${BCFTOOLS}')
OUT = f'{R}/l0'
os.makedirs(OUT, exist_ok=True)
byc = collections.defaultdict(set)
for line in open(f'{R}/deductive/l0_variants.tsv'):
    if line.startswith('gene\t'):
        continue
    k = line.split('\t')[1]
    byc[k.split(':')[0].replace('chr', '')].add(k)
tot = 0
for ch in sorted(byc, key=lambda x: int(x)):
    keys = byc[ch]
    if os.path.exists(f'{OUT}/chr{ch}.ds.done'):
        tot += len(keys)
        print(f'chr{ch}: skip (done)')
        continue
    with open(f'{OUT}/chr{ch}.reg', 'w') as fh:
        for k in sorted(keys, key=lambda x: int(x.split(':')[1])):
            p = k.split(':')[1]
            fh.write(f'{ch}\t{p}\t{p}\n')
    cmd = f"{BCF} query -R {OUT}/chr{ch}.reg -f '%CHROM:%POS:%REF:%ALT[\t%DS]\n' {R}/band_vcf/chr{ch}.band.vcf.gz > {OUT}/chr{ch}.ds.tsv"
    rc = subprocess.run(cmd, shell=True).returncode
    got = set()
    n = 0
    for line in open(f'{OUT}/chr{ch}.ds.tsv'):
        key = 'chr' + line.split('\t', 1)[0]
        n += 1
        if key in keys:
            got.add(key)
    miss = len(keys) - len(got)
    print(f'chr{ch}: want {len(keys):,} rows_fetched {n:,} matched {len(got):,} missing {miss}', flush=True)
    assert rc == 0 and miss == 0, f'GATE FAIL chr{ch}: rc={rc} missing={miss}'
    open(f'{OUT}/chr{ch}.ds.done', 'w').write('ok\n')
    tot += len(keys)
print(f'TOTAL variants extracted: {tot:,}')
print('L0_DS_COMPLETE')
