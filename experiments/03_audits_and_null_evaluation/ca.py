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
import csv, subprocess, glob, os, collections
B = _configured('${PHENO_DIR}')
BCF = _configured('${BCFTOOLS}')
gt = {s.split('_')[0] for s in subprocess.run(_configured(f'{BCF} query -l ${{PROJECT_ROOT}}/work/ref/band_vcf/chr22.band.vcf.gz'), shell=True, capture_output=True, text=True).stdout.split()}

def hdr(p):
    with open(p, encoding='utf-8', errors='ignore') as fh:
        return fh.readline().rstrip('\n').split('\t')
for rel in ['03.CA/CA_09_biochemical_analysis.tsv', '03.CA/CA_03_Medical_History.tsv', '03.CA/CA_01_subject_demographics.tsv']:
    h = hdr(_source_path(B, rel))
    keep = [c2 for c2 in h if any((k in c2.upper() for k in ('TCHL', 'CHOL', 'HDL', 'LDL', 'TG', 'GLU', 'AGE', 'SEX', 'HTN', 'DM', 'LIP', 'HYPER', 'DIAB')))]
    print(f'[{rel}] {len(h)} cols | relevant: {keep[:14]}')
print()
prs = sorted(glob.glob(f'{B}/other/PR/*.tsv'))
print('PR tables:', [os.path.basename(p) for p in prs][:10])
for p in prs[:3]:
    h = hdr(p)
    keep = [c2 for c2 in h if any((k in c2.upper() for k in ('TCHL', 'CHOL', 'AGE', 'SEX', 'HTN', 'DM')))]
    print(f'  [{os.path.basename(p)}] {len(h)} cols | relevant: {keep[:10]}')
print()
oth = sorted({os.path.dirname(os.path.relpath(p, B)) for p in glob.glob(f'{B}/other/**/*.tsv', recursive=True)})
print('other/ subdirs:', oth)
