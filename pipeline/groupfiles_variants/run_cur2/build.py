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
EIG = required('EIGENVEC_FILE')
OUT = os.environ.get('PHENO_CUR_OUT') or os.path.join(required('PROJECT_ROOT'), 'work/ref/pheno_cur')
os.makedirs(OUT, exist_ok=True)
pcs = {}
with open(EIG) as fh:
    fh.readline()
    for line in fh:
        fields = line.split()
        if len(fields) < 12:
            raise ValueError('Each eigenvector row must contain IID and ten PCs')
        pcs[fields[1]] = fields[2:12]
assert len(pcs) == N_SAMPLES, len(pcs)
missing = json_values('MISSING_VALUES_JSON')
sex, age = ({}, {})
for cfg in configs:
    prefix = cfg['prefix']
    subject = source_path(required(prefix + 'SUBJECT_FILE'))
    sex_column, age_column = (required(prefix + 'SEX_COLUMN'), required(prefix + 'AGE_COLUMN'))
    male, female = (json_values(prefix + 'MALE_VALUES_JSON'), json_values(prefix + 'FEMALE_VALUES_JSON'))
    if male & female:
        raise ValueError(f"Male/female encodings overlap for cohort {cfg['scope']}")
    with open(subject, newline='', errors='replace') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        if not {cfg['id_column'], sex_column, age_column}.issubset(reader.fieldnames or []):
            raise ValueError(f'Configured ID/sex/age columns are missing from {subject}')
        for row in reader:
            person = row[cfg['id_column']].strip()
            if person not in dup_of:
                continue
            sx, ag = ((row.get(sex_column) or '').strip(), (row.get(age_column) or '').strip())
            if sx in male:
                sex.setdefault(person, 1)
            elif sx in female:
                sex.setdefault(person, 0)
            if ag not in missing:
                age.setdefault(person, ag)
    print(f"cohort {cfg['scope']}: subject table loaded", flush=True)
TRAITS = ('HTN', 'DM', 'LIP', 'FRAC1', 'GASTRO', 'ARTH', 'ALLER', 'PER', 'FLIV', 'POL', 'OSTE', 'GASULCER', 'THY', 'CATA', 'MI', 'GB', 'BPH', 'ASTH', 'DUOULCER')
summary = {}
for trait in TRAITS:
    path = os.path.join(OUT, f'{trait.lower()}_cur.tsv')
    n_case = n_ctrl = 0
    with open(path, 'w') as fh:
        fh.write('sample_id\ty\tage\tsex_male\t' + '\t'.join((f'PC{i}' for i in range(1, 11))) + '\n')
        for person, genotype in dup_of.items():
            if genotype not in pcs or person not in sex or person not in age:
                continue
            state = status.get(trait, {}).get(person)
            if state is None:
                continue
            y = int(state == 'case')
            fh.write(f'{genotype}\t{y}\t{age[person]}\t{sex[person]}\t' + '\t'.join(pcs[genotype]) + '\n')
            n_case += y
            n_ctrl += 1 - y
    n = n_case + n_ctrl
    summary[trait] = {'cases': n_case, 'controls': n_ctrl, 'n': n, 'prevalence': round(n_case / n, 5) if n else None}
    print(f'{trait}: cases {n_case:,} controls {n_ctrl:,} n {n:,}', flush=True)
expected = json.loads(required('EXPECTED_CASE_COUNTS_JSON'))
if not isinstance(expected, dict) or not expected:
    raise ValueError('EXPECTED_CASE_COUNTS_JSON must be a nonempty mapping of trait names to case counts')
fails = []
for trait, count in expected.items():
    if not isinstance(count, int) or isinstance(count, bool) or count < 0:
        raise ValueError('Expected case counts must be nonnegative integers')
    got = summary.get(trait, {}).get('cases')
    if got != count:
        fails.append(f'{trait}: observed {got}, expected {count}')
with open(os.path.join(OUT, 'pheno_cur_summary.json'), 'w') as fh:
    json.dump({'summary': summary, 'expected_case_counts': expected, 'gate_failures': fails, 'note': 'Source columns and encodings supplied in configuration'}, fh, indent=1)
if fails:
    print('GATE FAIL:', fails)
    raise SystemExit(9)
print('GATE OK')
