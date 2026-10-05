#!/usr/bin/env python3

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
import os, glob, json, csv, gzip, subprocess, collections
B = _configured('${PHENO_DIR}')
OUT = _configured('${PROJECT_ROOT}/work/run_audit')
os.makedirs(OUT, exist_ok=True)
rep = {}
gt = set()
sam = _configured('${PROJECT_ROOT}/work/ref/band_vcf/chr22.band.vcf.gz')
BCF = _configured('${BCFTOOLS}')
for s in subprocess.run(f'{BCF} query -l {sam}', shell=True, capture_output=True, text=True).stdout.split():
    gt.add(s.strip())
print(f'genotyped samples in band VCF: {len(gt):,}', flush=True)
assert gt, 'GATE FAIL: no sample IDs read from band VCF'
exam = {}
for scope in 'ABC':
    exam[f'COHORT_{scope}'] = _source_path(B, _required(f'COHORT_{scope}_SUBJECT_FILE'))
print('subject tables:', {k: os.path.relpath(v, B) for k, v in exam.items()}, flush=True)
cohorts = {}
for coh, p in sorted(exam.items()):
    ids = set()
    with open(p, encoding='utf-8', errors='ignore') as fh:
        rd = csv.reader(fh, delimiter='\t')
        hdr = next(rd)
        source_id_index = _id_index(fh.name, hdr)
        idcol = _id_index(fh.name, hdr)
        for row in rd:
            if row:
                ids.add(row[idcol].strip())
    cohorts[coh] = {'table': os.path.relpath(p, B), 'n_rows': len(ids), 'n_genotyped': len(ids & gt), 'header_id_col': hdr[idcol]}
    print(f'  {coh:<5} rows {len(ids):>7,} | genotyped {len(ids & gt):>7,}', flush=True)
rep['Q1_cohorts'] = cohorts
rep['Q1_total_genotyped'] = len(gt)
rep['Q1_genotyped_covered'] = sum((c['n_genotyped'] for c in cohorts.values()))
TRAITS = json.loads(_required('TRAIT_COLUMN_KEYWORDS_JSON'))
found = collections.defaultdict(dict)
for p in glob.glob(f'{B}/**/*.tsv', recursive=True):
    base = os.path.basename(p)
    coh = base.split('_')[0]
    try:
        with open(p, encoding='utf-8', errors='ignore') as fh:
            hdr = fh.readline().rstrip('\n').split('\t')
    except OSError:
        continue
    up = [h.strip() for h in hdr]
    for tname, keys in TRAITS.items():
        for h in up:
            hu = h.upper()
            for k in keys:
                if hu.endswith('_' + k.upper()) or hu == k.upper() or hu.endswith(k.upper()):
                    found[tname].setdefault(coh, []).append(f'{base}:{h}')
                    break
for tname in TRAITS:
    per = {c: v[:3] for c, v in found[tname].items()}
    print(f'[{tname}] cohorts with a matching column: {sorted(per)}', flush=True)
rep['Q2_trait_columns'] = {t: {c: v for c, v in found[t].items()} for t in TRAITS}
cur = {}
cands = {'htn_v2': _configured('${PROJECT_ROOT}/work/ref/pheno_v2/htn_v2.tsv'), 'dm_v2': _configured('${PROJECT_ROOT}/work/ref/pheno_v2/dm_v2.tsv'), 'lip_v2': _configured('${PROJECT_ROOT}/work/ref/pheno_v2/lip_v2.tsv'), 'tchl_primary': _configured('${ANALYSIS_ROOT}/work/phenotype/tchl_primary_complete_case.tsv')}
coh_ids = {}
for coh, p in exam.items():
    with open(p, encoding='utf-8', errors='ignore') as fh:
        rd = csv.reader(fh, delimiter='\t')
        source_header = next(rd)
        source_id_index = _id_index(fh.name, source_header)
        coh_ids[coh] = {r[source_id_index].strip() for r in rd if r}
for name, p in cands.items():
    if not os.path.exists(p):
        cur[name] = {'exists': False}
        print(f'[{name}] MISSING', flush=True)
        continue
    with open(p, encoding='utf-8', errors='ignore') as fh:
        rd = csv.reader(fh, delimiter='\t')
        hdr = next(rd)
        source_id_index = _id_index(fh.name, hdr)
        ids = {r[source_id_index].strip() for r in rd if r}
    per = {c: len(ids & s) for c, s in coh_ids.items()}
    cur[name] = {'exists': True, 'n': len(ids), 'columns': hdr, 'by_cohort': per, 'n_genotyped': len(ids & gt)}
    print(f'[{name}] n={len(ids):,} genotyped={len(ids & gt):,} by_cohort={per}', flush=True)
    print(f"    columns: {','.join(hdr)}", flush=True)
rep['Q3_current_tables'] = cur
json.dump(rep, open(f'{OUT}/audit.json', 'w'), indent=1, ensure_ascii=False)
print('\nAUDIT_PART1_COMPLETE', flush=True)
