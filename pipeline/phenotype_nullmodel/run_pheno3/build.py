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
import math
TCHL_MIN = float(required("TCHL_MIN_VALUE"))
TCHL_MAX = float(required("TCHL_MAX_VALUE"))
if not math.isfinite(TCHL_MIN) or not math.isfinite(TCHL_MAX) or not TCHL_MIN < TCHL_MAX:
    raise ValueError("TCHL bounds must be finite and ordered")
import os, csv, glob, json, math, subprocess, collections, sys

BASE = required("PHENO_DIR")
EIG = required("EIGENVEC_FILE")
BCF = os.environ.get("BCFTOOLS", "bcftools")
VCF22 = required("GENOTYPE_FILE")
OUT = required("OUTPUT_DIR")
NPC  = 5
V2 = {t: int(required("PHENO_V2_" + t + "_CASES")) for t in ("HTN", "DM", "LIP")}
if any(value < 0 for value in V2.values()): raise ValueError("V2 case counts must be nonnegative")
MIN_TCHL_SAMPLES = int(required("MIN_TCHL_SAMPLES"))
if MIN_TCHL_SAMPLES <= 0: raise ValueError("MIN_TCHL_SAMPLES must be positive")
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
        f = line.split(); pcs[f[1]] = f[2:2+NPC]
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

BIN = ("HTN","DM","LIP")
dis = {t: {} for t in BIN}
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

TCH = [(source_path(required("COHORT_" + cfg["tag"] + "_LAB_FILE")),
        required("COHORT_" + cfg["tag"] + "_TCHL_COLUMN"),
        required("COHORT_" + cfg["tag"] + "_LAB_ID_COLUMN")) for cfg in cohorts]
tchl = {}
for path, col, id_col in TCH:
    d = read_col(path, col, id_col)
    assert d, f"GATE FAIL: {col} yielded nothing from {path}"
    n0 = len(tchl)
    for k, v in d.items():
        try: x = float(v)
        except ValueError: continue
        if TCHL_MIN < x < TCHL_MAX: tchl.setdefault(k, x)
    print(f"  {col}: {len(d):,} rows -> +{len(tchl)-n0:,} usable", flush=True)
print(f"TCHL pooled: {len(tchl):,}", flush=True)

def _erfinv(y):
    a = 0.147
    ln = math.log(1 - y*y) if abs(y) < 1 else -700.0
    t = 2/(math.pi*a) + ln/2
    return math.copysign(math.sqrt(max(math.sqrt(t*t - ln/a) - t, 0.0)), y)
def rint(d):
    items = sorted(d.items(), key=lambda kv: kv[1]); n = len(items); out = {}; i = 0
    while i < n:
        j = i
        while j+1 < n and items[j+1][1] == items[i][1]: j += 1
        r = (i+j)/2 + 1
        p = (r - 0.375)/(n + 0.25)
        z = math.sqrt(2)*_erfinv(2*p - 1)
        for k in range(i, j+1): out[items[k][0]] = z
        i = j+1
    return out
tchl_z = rint(tchl)

TRAITS = ["HTN","DM","LIP","TCHL"]
summary = {}
for t in TRAITS:
    path = f"{OUT}/{t.lower()}_v3.tsv"
    n_case = n = 0; byc = collections.Counter()
    with open(path,"w") as fh:
        fh.write("sample_id\ty\tage\tsex_male\tCT\tNC\t" +
                 "\t".join(f"PC{i}" for i in range(1,NPC+1)) + "\n")
        for nk, dk in dup_of.items():
            if dk not in pcs or nk not in sex or nk not in age or nk not in cohort: continue
            if t == "TCHL":
                if nk not in tchl_z: continue
                y = f"{tchl_z[nk]:.6f}"; is_case = 0
            else:
                st = dis[t].get(nk)
                if st is None: continue
                is_case = 1 if st == "case" else 0; y = str(is_case)
            sm = 1 if sex[nk] in MALE_VALUES else 0
            ct = 1 if cohort[nk] == "A" else 0
            nc = 1 if cohort[nk] == "B" else 0
            fh.write(f"{dk}\t{y}\t{age[nk]}\t{sm}\t{ct}\t{nc}\t" + "\t".join(pcs[dk]) + "\n")
            n += 1; n_case += is_case; byc[cohort[nk]] += 1
    summary[t] = {"n": n, "cases": (n_case if t != "TCHL" else None),
                  "by_cohort": dict(byc),
                  "trait_type": "quantitative" if t == "TCHL" else "binary"}
    print(f"[{t}] n={n:,} cases={n_case if t!='TCHL' else '-'} by_cohort={dict(byc)}", flush=True)

json.dump(summary, open(f"{OUT}/pheno_v3_summary.json","w"), indent=1)


fails = []
if len(summary) != 4: fails.append(f"wrote {len(summary)} traits, expected 4")
for t, e in V2.items():
    got = summary[t]["cases"]
    if got != e: fails.append(f"{t}: cases {got:,} != v2 {e:,} — binary definition drifted")
for t, s in summary.items():
    if s["n"] == 0: fails.append(f"{t}: zero rows")
    if set(s["by_cohort"]) != {"A","B","C"}:
        fails.append(f"{t}: cohorts {sorted(s['by_cohort'])} != A/B/C")
if summary["TCHL"]["n"] < MIN_TCHL_SAMPLES:
    fails.append(f"TCHL n={summary['TCHL']['n']:,} — below configured minimum over 3 cohorts")
if fails:
    for m in fails: print("GATE FAIL:", m)
    sys.exit(6)
print("PHENO_V3_COMPLETE")
