
def _id_index(path, header):
    path = _config_os.path.abspath(path)
    for scope in "ABCD":
        prefix = f"COHORT_{scope}_"
        for suffix in ("SUBJECT_FILE", "BIOCHEM_FILE", "DISEASE_FILE"):
            configured = _config_os.environ.get(prefix + suffix)
            if configured:
                base = _config_os.environ.get("PHENO_DIR", "")
                if _config_os.path.abspath(_source_path(base, configured)) == path:
                    return header.index(_required(prefix + "ID_COLUMN"))
    generic = _config_os.environ.get("ID_COLUMN")
    if generic and generic in header:
        return header.index(generic)
    for column in ("sample_id", "SAMPLE_ID", "IID", "ID"):
        if column in header:
            return header.index(column)
    raise ValueError("Set ID_COLUMN for source tables without a scoped configuration")
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
import csv, glob, os, subprocess, collections
B = _configured('${PHENO_DIR}')
BCF = _configured('${BCFTOOLS}')
GT = {s.strip() for s in subprocess.run(_configured(f'{BCF} query -l ${{PROJECT_ROOT}}/work/ref/band_vcf/chr22.band.vcf.gz'), shell=True, capture_output=True, text=True).stdout.split()}
bare = {g.split('_')[0] for g in GT}
PC = _configured('${EIGENVEC_FILE}')
print('PC file exists:', os.path.exists(PC))
if os.path.exists(PC):
    with open(PC) as fh:
        h = fh.readline().split()
    print('  PC header:', h[:14])
    n = sum((1 for _ in open(PC))) - 1
    print(f'  rows: {n:,}')
    with open(PC) as fh:
        next(fh)
        ids = {l.split()[1] if len(l.split()) > 1 else l.split()[0] for l in fh}
    print(f'  IDs matching genotypes (col2/IID): {len(ids & GT):,} of {len(GT):,}')
E = {}
for coh, rel in [('CT', _configured('${COHORT_A_SUBJECT_FILE}')), ('AS', _configured('${COHORT_C_SUBJECT_FILE}')), ('NC', _configured('${COHORT_B_SUBJECT_FILE}'))]:
    with open(_source_path(B, rel), encoding='utf-8', errors='ignore') as fh:
        rd = csv.reader(fh, delimiter='\t')
        source_header = next(rd)
        source_id_index = _id_index(fh.name, source_header)
        E[coh] = {r[source_id_index].strip() for r in rd if r}
print()
asg = collections.Counter()
for s in bare:
    hits = [c2 for c2 in E if s in E[c2]]
    asg['+'.join(hits) if hits else 'NONE'] += 1
print('cohort assignment of genotyped samples:')
for k, v in asg.most_common(8):
    print(f"  {k or 'NONE':<10} {v:,}")
print()
for coh, rel, pre in [('CT', _configured('${COHORT_A_SUBJECT_FILE}'), 'CT'), ('AS', _configured('${COHORT_C_SUBJECT_FILE}'), 'COHORT_C'), ('NC', _configured('${COHORT_B_SUBJECT_FILE}'), 'NC')]:
    with open(_source_path(B, rel), encoding='utf-8', errors='ignore') as fh:
        h = fh.readline().rstrip('\n').split('\t')
    ag = [c2 for c2 in h if 'AGE' in c2.upper()]
    sx = [c2 for c2 in h if 'SEX' in c2.upper()]
    print(f'  {coh}: age={ag[:3]} sex={sx[:3]}')
print()
for coh, rel, col in [('CT', _configured('${COHORT_A_BIOCHEM_FILE}'), _required('COHORT_A_TCHL_COLUMN')), ('NC', _configured('${COHORT_B_BIOCHEM_FILE}'), _required('COHORT_B_TCHL_COLUMN')), ('AS', _configured('${COHORT_C_BIOCHEM_FILE}'), _required('COHORT_C_TCHL_COLUMN'))]:
    with open(_source_path(B, rel), encoding='utf-8', errors='ignore') as fh:
        h = fh.readline().rstrip('\n').split('\t')
    print(f"  {coh}: '{col}' present={col in h}")
