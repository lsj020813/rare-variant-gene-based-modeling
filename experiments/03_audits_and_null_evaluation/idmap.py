
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
gt = {s.strip() for s in subprocess.run(_configured(f'{BCF} query -l ${{PROJECT_ROOT}}/work/ref/band_vcf/chr22.band.vcf.gz'), shell=True, capture_output=True, text=True).stdout.split()}

def bare(x):
    return x.split('_')[0]
gt_bare = {bare(s) for s in gt}
print(f'genotyped: {len(gt):,} | bare-unique: {len(gt_bare):,}')
print('genotype ID prefixes:', collections.Counter((s[:6] for s in sorted(gt_bare))).most_common(6))
exam = {}
for scope in 'ABC':
    exam[f'COHORT_{scope}'] = _source_path(B, _required(f'COHORT_{scope}_SUBJECT_FILE'))
print()
tot = 0
for coh, p in sorted(exam.items()):
    with open(p, encoding='utf-8', errors='ignore') as fh:
        rd = csv.reader(fh, delimiter='\t')
        source_header = next(rd)
        source_id_index = _id_index(fh.name, source_header)
        ids = {r[source_id_index].strip() for r in rd if r}
    ov = len(ids & gt_bare)
    tot += ov
    pre = collections.Counter((i[:6] for i in list(ids)[:2000])).most_common(2)
    print(f'  {coh:<5} rows {len(ids):>7,} | genotyped {ov:>7,} | id prefixes {pre}')
print(f'\nsum of per-cohort genotyped: {tot:,} vs total genotyped {len(gt_bare):,}')
for name, p in [('htn_v2', _configured('${PROJECT_ROOT}/work/ref/pheno_v2/htn_v2.tsv')), ('tchl_primary', _configured('${ANALYSIS_ROOT}/work/phenotype/tchl_primary_complete_case.tsv'))]:
    with open(p, encoding='utf-8', errors='ignore') as fh:
        rd = csv.reader(fh, delimiter='\t')
        source_header = next(rd)
        source_id_index = _id_index(fh.name, source_header)
        ids = {bare(r[source_id_index]) for r in rd if r}
    print(f'\n[{name}] n={len(ids):,} in-genotype={len(ids & gt_bare):,}')
    for coh, pp in sorted(exam.items()):
        with open(pp, encoding='utf-8', errors='ignore') as fh:
            rd = csv.reader(fh, delimiter='\t')
            source_header = next(rd)
            source_id_index = _id_index(fh.name, source_header)
            cids = {r[source_id_index].strip() for r in rd if r}
        print(f'    {coh}: {len(ids & cids):,}')
