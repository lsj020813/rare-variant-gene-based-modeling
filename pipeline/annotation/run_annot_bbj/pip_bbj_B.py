#!/usr/bin/env python3
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import os, sys, gzip, json, argparse, collections
L1 = _config_path("${PROJECT_ROOT}/work/ref/annot/bbj_l1")
ap = argparse.ArgumentParser(); ap.add_argument("--annot", default="bbj_annot.tsv.gz"); ap.add_argument("--outdir", default="bbj_pip"); ap.add_argument("--traits-out", default="traits.tsv"); ap.add_argument("--fold-map", default=None, help="JSON {chr: fold} replacing the default v10 FOLDS (subset runs with --folds-override)")
args = ap.parse_args()
GROUP = {"TC": "지질", "LDLC": "지질", "HDLC": "지질", "TG": "지질",
         "T2D": "당대사", "HbA1c": "당대사", "Glucose": "당대사",
         "SBP": "혈압", "DBP": "혈압", "MAP": "혈압", "PP": "혈압",
         "RBC": "혈액", "Hb": "혈액", "Ht": "혈액", "MCV": "혈액", "MCH": "혈액", "MCHC": "혈액", "WBC": "혈액", "NEU": "혈액", "LYM": "혈액", "MON": "혈액", "EOS": "혈액", "BAS": "혈액", "PLT": "혈액"}
FOLD = {7:0,13:0,16:0,19:0,21:0, 2:1,8:1,9:1,10:1,15:1, 1:2,6:2,14:2, 3:3,12:3,17:3,22:3, 4:4,5:4,11:4,18:4,20:4}
if args.fold_map: FOLD = {int(k): int(v) for k, v in json.loads(args.fold_map).items()}; print(f"[B] fold map override: {FOLD}")
keys = set(); chrs_present = set()
with gzip.open(f"{L1}/{args.annot}", "rt") as fh:
    next(fh)
    for line in fh:
        k = line.split("\t", 1)[0]; keys.add(k); chrs_present.add(int(k.split(":")[0]))
folds_present = {FOLD[c] for c in chrs_present}
print(f"[B] annot keys {len(keys):,} chrs {sorted(chrs_present)} folds {sorted(folds_present)}", flush=True)
od = f"{L1}/{args.outdir}"; os.makedirs(od, exist_ok=True)
raw = sorted(f[:-4] for f in os.listdir(f"{L1}/raw/pip_raw") if f.endswith(".tsv"))
assert len(raw) == 76, f"GATE FAIL: traits {len(raw)}"
rows = []; tot_kept = 0; tot_all = 0
for tr in raw:
    seen = set(); kept = n = 0; regions = set(); cs = set(); csn = 0; pipmass = 0.0; folds_seen = set()
    with open(f"{L1}/raw/pip_raw/{tr}.tsv") as fh, gzip.open(f"{od}/{tr}.tsv.gz.tmp", "wt") as o:
        o.write(next(fh))
        for line in fh:
            n += 1; k = line.split("\t", 1)[0]
            if k not in keys: continue
            assert k not in seen, f"GATE FAIL: duplicate (trait,variant) {tr}"
            seen.add(k); kept += 1; o.write(line)
            a = line.rstrip("\n").split("\t"); regions.add(a[7]); pipmass += float(a[5]); folds_seen.add(FOLD[int(a[1])])
            if a[6] != "-1": cs.add((a[7], a[6])); csn += 1
    os.replace(f"{od}/{tr}.tsv.gz.tmp", f"{od}/{tr}.tsv.gz")
    tot_kept += kept; tot_all += n
    rows.append(dict(trait=tr, trait_group=GROUP.get(tr, "기타"), N="NA", rows_raw=n, rows_kept=kept, kept_pct=round(kept / n * 100, 2) if n else 0, regions=len(regions), cs=len(cs), cs_rows=csn, pip_mass=round(pipmass, 2), folds_seen_set=folds_seen, group_star="" if tr in GROUP else "★"))
    print(f"[B] {tr}: {kept:,}/{n:,} 행 ({rows[-1]['kept_pct']}%) 구역 {len(regions)} CS {len(cs)} CS행 {csn:,} ΣPIP {pipmass:.1f}", flush=True)
empty = [r["trait"] for r in rows if not (r["folds_seen_set"] >= folds_present)]
for r in rows:
    r["folds_seen"] = "".join(str(f) for f in sorted(r["folds_seen_set"])); del r["folds_seen_set"]
    if r["trait"] in empty: os.remove(f"{od}/{r['trait']}.tsv.gz")
with open(f"{L1}/{args.traits_out}.tmp", "w") as o:
    o.write("trait\ttrait_group\tN\n")
    for r in rows:
        if r["trait"] not in empty: o.write(f"{r['trait']}\t{r['trait_group']}\t{r['N']}\n")
print(f"[B] traits excluded (missing some chromosome fold): {len(empty)} {empty}")
os.replace(f"{L1}/{args.traits_out}.tmp", f"{L1}/{args.traits_out}")
with open(f"{L1}/{args.traits_out}.audit.tsv", "w") as o:
    o.write("\t".join(rows[0].keys()) + "\n")
    for r in rows: o.write("\t".join(str(v) for v in r.values()) + "\n")
grp = collections.Counter(r["trait_group"] for r in rows if r["trait"] not in empty)
open(f"{od}/.done", "w").write(f"ok {tot_kept}\n")
print(json.dumps(dict(traits=len(rows), traits_kept=len(rows)-len(empty), rows_kept=tot_kept, rows_raw=tot_all, groups=grp), ensure_ascii=False)); print("PIP_B_DONE")
