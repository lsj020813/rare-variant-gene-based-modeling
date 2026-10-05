
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
import csv, subprocess, collections, statistics as st
B = _configured('${PHENO_DIR}')
BCF = _configured('${BCFTOOLS}')
gt = {s.split('_')[0] for s in subprocess.run(_configured(f'{BCF} query -l ${{PROJECT_ROOT}}/work/ref/band_vcf/chr22.band.vcf.gz'), shell=True, capture_output=True, text=True).stdout.split()}

def col_stats(rel, col):
    p = _source_path(B, rel)
    with open(p, encoding='utf-8', errors='ignore') as fh:
        rd = csv.reader(fh, delimiter='\t')
        hdr = next(rd)
        source_id_index = _id_index(fh.name, hdr)
        if col not in hdr:
            return None
        ci = hdr.index(col)
        vals = []
        big = collections.Counter()
        ng = 0
        for r in rd:
            if not r or len(r) <= ci:
                continue
            if r[source_id_index].strip() not in gt:
                continue
            ng += 1
            v = r[ci].strip()
            if v in ('', 'NA', '.'):
                continue
            try:
                f = float(v)
            except ValueError:
                continue
            if f >= 600 or f <= 0:
                big[v] += 1
            else:
                vals.append(f)
        return dict(genotyped=ng, valid=len(vals), outliers=big.most_common(5), mean=round(st.mean(vals), 1) if vals else None, sd=round(st.pstdev(vals), 1) if vals else None, lo=min(vals) if vals else None, hi=max(vals) if vals else None)
for rel, col in [(_configured('${COHORT_A_BIOCHEM_FILE}'), _required('COHORT_A_TCHL_COLUMN')), (_configured('${COHORT_B_BIOCHEM_FILE}'), _required('COHORT_B_TCHL_COLUMN')), (_configured('${COHORT_C_BIOCHEM_FILE}'), _required('COHORT_C_TCHL_ORIGINAL_COLUMN')), (_configured('${COHORT_C_BIOCHEM_FILE}'), _required('COHORT_C_TCHL_COLUMN'))]:
    s = col_stats(rel, col)
    print(f'[{col}] {s}')
