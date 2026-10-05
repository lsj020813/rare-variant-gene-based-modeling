
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
import csv, subprocess, glob, os, collections
B = _configured('${PHENO_DIR}')
BCF = _configured('${BCFTOOLS}')
gt = {s.split('_')[0] for s in subprocess.run(_configured(f'{BCF} query -l ${{PROJECT_ROOT}}/work/ref/band_vcf/chr22.band.vcf.gz'), shell=True, capture_output=True, text=True).stdout.split()}

def ids(p):
    try:
        with open(p, encoding='utf-8', errors='ignore') as fh:
            rd = csv.reader(fh, delimiter='\t')
            source_header = next(rd)
            source_id_index = _id_index(fh.name, source_header)
            return {r[source_id_index].strip() for r in rd if r}
    except Exception:
        return set()
covered = set()
for p in [_required('COHORT_A_SUBJECT_FILE'), _required('COHORT_B_SUBJECT_FILE'), _required('COHORT_C_SUBJECT_FILE')]:
    covered |= ids(p)
missing = gt - covered
print(f'unmapped genotyped: {len(missing):,}')
best = []
for p in sorted(glob.glob(f'{B}/**/*.tsv', recursive=True)):
    s = ids(p)
    ov = len(s & missing)
    if ov > 100:
        best.append((ov, os.path.relpath(p, B), len(s)))
best.sort(reverse=True)
print('tables containing the unmapped samples (top 12):')
for ov, rel, n in best[:12]:
    print(f'  {ov:>7,} / {n:>7,}  {rel}')
if not best:
    print('  NONE — these IDs appear in no phenotype table at all')
