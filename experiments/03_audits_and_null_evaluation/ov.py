
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
import csv, subprocess, glob, os
B = _configured('${PHENO_DIR}')
BCF = _configured('${BCFTOOLS}')
gt = {s.split('_')[0] for s in subprocess.run(_configured(f'{BCF} query -l ${{PROJECT_ROOT}}/work/ref/band_vcf/chr22.band.vcf.gz'), shell=True, capture_output=True, text=True).stdout.split()}

def ids(p):
    with open(p, encoding='utf-8', errors='ignore') as fh:
        rd = csv.reader(fh, delimiter='\t')
        source_header = next(rd)
        source_id_index = _id_index(fh.name, source_header)
        return {r[source_id_index].strip() for r in rd if r}
E = {c2: ids(p) for c2, p in [('COHORT_C', _source_path(B, _configured('${COHORT_C_SUBJECT_FILE}'))), ('COHORT_D', _source_path(B, _configured('${COHORT_D_SUBJECT_FILE}'))), ('CT', _source_path(B, _configured('${COHORT_A_SUBJECT_FILE}'))), ('NC', _source_path(B, _configured('${COHORT_B_SUBJECT_FILE}')))]}
print('pairwise overlaps (genotyped only):')
ks = sorted(E)
for i, a in enumerate(ks):
    for b2 in ks[i + 1:]:
        ov = len(E[a] & E[b2] & gt)
        if ov:
            print(f'  {a} n {b2}: {ov:,}')
print()
print('COHORT_D subset of COHORT_C?', E['COHORT_D'] & gt <= E['COHORT_C'] & gt)
print('COHORT_C genotyped:', f"{len(E['COHORT_C'] & gt):,}", '| COHORT_D genotyped:', f"{len(E['COHORT_D'] & gt):,}")
union = set().union(*[E[k] & gt for k in ks])
print('union of 4 cohorts (genotyped):', f'{len(union):,}', '| genotyped total:', f'{len(gt):,}')
print('genotyped NOT in any subject table:', f'{len(gt - union):,}')
rest = sorted(gt - union)[:5]
print('examples:', rest)
import collections
print('their prefixes:', collections.Counter((x[:6] for x in gt - union)).most_common(5))
print()
print('COHORT_D tables:', sorted((os.path.basename(p) for p in glob.glob(_configured('${COHORT_D_TABLE_GLOB}'))))[:8])
