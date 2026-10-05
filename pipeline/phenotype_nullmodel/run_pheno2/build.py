#!/usr/bin/env python3
import os
import json

def required(name):
    value = os.environ.get(name, "")
    if not value.strip():
        raise RuntimeError("Required environment variable is missing or empty: " + name)
    return value

N_SAMPLES = int(required("N_SAMPLES"))
if N_SAMPLES <= 0: raise ValueError("N_SAMPLES must be positive")
CASE_VALUES = set(json.loads(required("CASE_VALUES_JSON")))
CONTROL_VALUES = set(json.loads(required("CONTROL_VALUES_JSON")))
MALE_VALUES = set(json.loads(required("MALE_VALUES_JSON")))
if not CASE_VALUES or not CONTROL_VALUES or CASE_VALUES & CONTROL_VALUES:
    raise ValueError("Case and control value sets must be nonempty and disjoint")
if not MALE_VALUES or not all(isinstance(v, str) for v in CASE_VALUES | CONTROL_VALUES | MALE_VALUES):
    raise ValueError("Declared phenotype codes must be nonempty sets of strings")
import os, glob, json, subprocess, collections

BCF = os.environ.get("BCFTOOLS", "bcftools")
VCF22 = required("GENOTYPE_FILE")
BASE = required("PHENO_DIR")
EIG = required("EIGENVEC_FILE")
OUT = required("OUTPUT_DIR")
V1 = {t: int(required("PHENO_V1_" + t + "_CASES")) for t in ("HTN", "DM", "LIP")}
if any(value < 0 for value in V1.values()): raise ValueError("V1 case counts must be nonnegative")
ADDED_MIN = {t: int(required("PHENO_ADDED_COHORT_MIN_" + t + "_CASES")) for t in ("HTN", "DM", "LIP")}
if any(value < 0 for value in ADDED_MIN.values()): raise ValueError("ADDED_MIN case counts must be nonnegative")
os.makedirs(OUT, exist_ok=True)

raw = subprocess.run([BCF,"query","-l",VCF22],capture_output=True,text=True,check=True).stdout.split()
def norm(s):
    h = s.split("_"); return h[0] if len(h)==2 and h[0]==h[1] else s
dup_of = {norm(s): s for s in raw}
assert len(dup_of) == N_SAMPLES, len(dup_of)

pcs = {}
with open(EIG) as fh:
    fh.readline()
    for line in fh:
        f = line.split(); pcs[f[1]] = f[2:12]
assert len(pcs) == N_SAMPLES, len(pcs)

MISS = set(json.loads(required("MISSING_CODES_JSON")))
def read_col(path, col, id_col):
    d = {}
    if not os.path.isfile(path): raise FileNotFoundError(path)
    with open(path, newline="", errors="replace") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        if col not in hdr or id_col not in hdr: raise ValueError("Configured column missing in " + path)
        i = hdr.index(col); iid = hdr.index(id_col)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) != len(hdr): continue
            k, v = p[iid], p[i].strip()
            if k in dup_of and v not in MISS: d.setdefault(k, v)
    return d

cohorts = []
for tag in ("A", "B", "C"):
    prefix = "COHORT_" + tag + "_"
    cohort_config = dict(tag=tag,
        subject_file=required(prefix + "SUBJECT_FILE"),
        subject_id=required(prefix + "SUBJECT_ID_COLUMN"),
        sex_column=required(prefix + "SEX_COLUMN"),
        age_column=required(prefix + "AGE_COLUMN"),
        disease_glob=required(prefix + "DISEASE_GLOB"),
        disease_id=required(prefix + "DISEASE_ID_COLUMN"),
        disease_columns={trait: json.loads(required(prefix + trait + "_COLUMNS_JSON")) for trait in ("HTN", "DM", "LIP")})
    if any(not isinstance(v, list) or not v or not all(isinstance(c, str) and c for c in v) for v in cohort_config["disease_columns"].values()):
        raise ValueError("Diagnosis column lists must contain declared source columns for cohort " + tag)
    cohorts.append(cohort_config)

def source_path(value):
    return value if os.path.isabs(value) else os.path.join(BASE, value)

sex, age, cohort = {}, {}, {}
for cfg in cohorts:
    path = source_path(cfg["subject_file"])
    s = read_col(path, cfg["sex_column"], cfg["subject_id"])
    a = read_col(path, cfg["age_column"], cfg["subject_id"])
    for k, v in s.items(): sex.setdefault(k, v)
    for k, v in a.items(): age.setdefault(k, v)
    for k in set(s) | set(a): cohort.setdefault(k, cfg["tag"])
    print(f"subject cohort {cfg['tag']}: sex {len(s):,} age {len(a):,}", flush=True)
cc = collections.Counter(cohort.values())
assert set(cc) == {"A", "B", "C"}, f"GATE FAIL: cohorts {sorted(cc)}"

TRAITS = ("HTN","DM","LIP")
dis = {t: {} for t in TRAITS}

def apply(trait, person, value):
    if value in CASE_VALUES: dis[trait][person] = "case"
    elif value in CONTROL_VALUES and dis[trait].get(person) != "case": dis[trait][person] = "ctrl"

for cfg in cohorts:
    files = [p for p in sorted(glob.glob(source_path(cfg["disease_glob"]), recursive=True)) if "donotuse" not in p]
    if not files: raise FileNotFoundError("No diagnosis tables for cohort " + cfg["tag"])
    print(f"disease tables cohort {cfg['tag']}: {len(files)}", flush=True)
    found_columns = set()
    for path in files:
        with open(path, newline="", errors="replace") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            if cfg["disease_id"] not in hdr: raise ValueError("Configured diagnosis ID column missing in " + path)
            iid = hdr.index(cfg["disease_id"])
            cols = {trait: [hdr.index(column) for column in names if column in hdr]
                    for trait, names in cfg["disease_columns"].items()}
            found_columns.update(column for names in cfg["disease_columns"].values() for column in names if column in hdr)
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) != len(hdr): continue
                k = p[iid]
                if k not in dup_of: continue
                for trait, indices in cols.items():
                    for i in indices: apply(trait, k, p[i].strip())
    expected_columns = {column for names in cfg["disease_columns"].values() for column in names}
    if found_columns != expected_columns:
        raise ValueError("Declared diagnosis columns missing for cohort " + cfg["tag"])

summary = {}
for t in TRAITS:
    path = f"{OUT}/{t.lower()}_v2.tsv"
    n_case = n_ctrl = 0
    with open(path,"w") as fh:
        fh.write("sample_id\ty\tage\tsex_male\t" + "\t".join(f"PC{i}" for i in range(1,11)) + "\n")
        for nk, dk in dup_of.items():
            if dk not in pcs or nk not in sex or nk not in age: continue
            st = dis[t].get(nk)
            if st is None: continue
            y = 1 if st == "case" else 0
            sex_male = 1 if sex[nk] in MALE_VALUES else 0
            fh.write(f"{dk}\t{y}\t{age[nk]}\t{sex_male}\t" + "\t".join(pcs[dk]) + "\n")
            n_case += y; n_ctrl += 1-y
    summary[t] = {"cases": n_case, "controls": n_ctrl, "n": n_case+n_ctrl,
                  "prevalence": round(n_case/(n_case+n_ctrl), 5)}
    print(f"{t}: cases {n_case:,} controls {n_ctrl:,} n {n_case+n_ctrl:,}", flush=True)



fails = []
for t, exp in V1.items():
    got = summary[t]["cases"]
    if got < exp + ADDED_MIN[t]:
        fails.append(f"{t}: {got} < {exp}+{ADDED_MIN[t]} — additional cohort not merged?")
json.dump({"summary": summary, "v1_reference": V1, "gate_failures": fails,
           "note": "v2 merges configured source cohorts with explicit diagnosis column mappings"},
          open(f"{OUT}/pheno_v2_summary.json","w"), indent=1)
if fails:
    print("GATE FAIL:", fails); raise SystemExit(9)
print("GATE OK — additional cohort merged")
