#!/usr/bin/env python3
import os
import glob
import json
import csv
import subprocess
import collections

def required(name):
    value = os.environ.get(name)
    if value is None or not value.strip():
        raise RuntimeError(f'Set {name} before running this script')
    return value

def required_int(name):
    value = int(required(name))
    if value <= 0:
        raise ValueError(f'{name} must be a positive integer')
    return value

def source_path(value):
    return value if os.path.isabs(value) else os.path.join(required('PHENO_DIR'), value)

def json_values(name):
    values = json.loads(required(name))
    if not isinstance(values, list) or not values:
        raise ValueError(f'{name} must be a nonempty JSON list')
    return {str(value) for value in values}

def cohort_config(scope):
    prefix = f'COHORT_{scope}_'
    mapping = json.loads(required(prefix + 'DISEASE_COLUMN_MAP_JSON'))
    if not isinstance(mapping, dict) or not mapping or (not all((isinstance(k, str) and isinstance(v, str) for k, v in mapping.items()))):
        raise ValueError(f'{prefix}DISEASE_COLUMN_MAP_JSON must map source columns to trait names')
    pattern = os.environ.get(prefix + 'DISEASE_GLOB')
    paths = sorted(glob.glob(source_path(pattern), recursive=True)) if pattern else [source_path(required(prefix + 'DISEASE_FILE'))]
    if not paths:
        raise ValueError(f'No disease files configured for cohort {scope}')
    cases = json_values(prefix + 'CASE_VALUES_JSON')
    controls = json_values(prefix + 'CONTROL_VALUES_JSON')
    if cases & controls:
        raise ValueError(f'Case/control encodings overlap for cohort {scope}')
    return {'scope': scope, 'prefix': prefix, 'paths': paths, 'id_column': required(prefix + 'ID_COLUMN'), 'mapping': mapping, 'cases': cases, 'controls': controls}

def norm(sample):
    fields = sample.split('_')
    return fields[0] if len(fields) == 2 and fields[0] == fields[1] else sample

def disease_status(configs, genotyped):
    status = {trait: {} for cfg in configs for trait in cfg['mapping'].values()}
    columns = collections.Counter()
    for cfg in configs:
        for path in cfg['paths']:
            with open(path, newline='', errors='replace') as fh:
                reader = csv.DictReader(fh, delimiter='\t')
                header = reader.fieldnames or []
                if cfg['id_column'] not in header:
                    raise ValueError(f'Missing configured ID column in {path}')
                selected = {column: trait for column, trait in cfg['mapping'].items() if column in header}
                if not selected:
                    raise ValueError(f'No configured disease columns found in {path}')
                columns.update(selected.values())
                for row in reader:
                    person = row[cfg['id_column']].strip()
                    if person not in genotyped:
                        continue
                    for column, trait in selected.items():
                        value = (row.get(column) or '').strip()
                        if value in cfg['cases']:
                            status[trait][person] = 'case'
                        elif value in cfg['controls'] and status[trait].get(person) != 'case':
                            status[trait][person] = 'ctrl'
    return (status, columns)
BCF = os.environ.get('BCFTOOLS', 'bcftools')
VCF22 = required('GENOTYPE_VCF_FILE')
N_SAMPLES = required_int('N_SAMPLES')
raw = subprocess.run([BCF, 'query', '-l', VCF22], capture_output=True, text=True, check=True).stdout.split()
dup_of = {norm(sample): sample for sample in raw}
if len(dup_of) != len(raw):
    raise ValueError('Genotype IDs collide after normalization')
assert len(dup_of) == N_SAMPLES, len(dup_of)
configs = [cohort_config(scope) for scope in 'ABC']
status, column_counts = disease_status(configs, dup_of)
OUT = os.environ.get('PHENO_CENSUS_OUT') or os.path.join(required('PROJECT_ROOT'), 'work/ref/pheno_census')
os.makedirs(OUT, exist_ok=True)
EVAL = ['HTN', 'DM', 'LIP']
eval_cases = {trait: {person for person, state in status.get(trait, {}).items() if state == 'case'} for trait in EVAL if trait in status}
rows = []
for trait, values in status.items():
    cases = {person for person, state in values.items() if state == 'case'}
    nc, nk = (len(cases), sum((state == 'ctrl' for state in values.values())))
    if not nc:
        continue
    overlap = {name: round(len(cases & persons) / nc, 4) for name, persons in eval_cases.items()}
    rows.append({'trait': trait, 'cases': nc, 'controls': nk, 'n_phenotyped': nc + nk, 'prevalence': round(nc / (nc + nk), 5), 'n_eff': round(4 * nc * nk / (nc + nk)), 'overlap_with_eval': overlap, 'n_columns': column_counts[trait]})
rows.sort(key=lambda row: -row['cases'])
with open(os.path.join(OUT, 'trait_census.json'), 'w') as fh:
    json.dump({'n_geno': len(dup_of), 'disease_files': sum((len(cfg['paths']) for cfg in configs)), 'definition': 'Configured case encoding takes precedence over configured control encoding across source tables; genotyped subset only', 'traits': rows}, fh, indent=1)
print(f'traits with >=1 case: {len(rows)}')
for row in rows[:45]:
    print(f"{row['trait']:<14} case {row['cases']:>7,} ctrl {row['controls']:>7,} prev {row['prevalence']} n_eff {row['n_eff']:>7,} cols {row['n_columns']}")
