#!/usr/bin/env python
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import gzip, os, json, subprocess, collections
P = _config_path("${PROJECT_ROOT}/work/fset/primary"); TABIX = "tabix"
SRC = _config_path("${PROJECT_ROOT}/work/ref/gpnmsa/scores.tsv.bgz")
RC = {"A": "T", "C": "G", "G": "C", "T": "A"}
def rc(s): return "".join(RC.get(b, "N") for b in reversed(s))
rows = []
with gzip.open(P + "/gpn_scores.tsv.gz", "rt") as f:
    hdr = f.readline()
    for l in f: rows.append(l.rstrip("\n").split("\t"))
mmpos = sorted({(r[1], int(r[2])) for r in rows if r[6] == "ref_mismatch"}, key=lambda x: (int(x[0]), x[1]))
bed = P + "/gpn_mismatch_positions.bed"
open(bed, "w").write("".join(f"{c}\t{p-1}\t{p}\n" for c, p in mmpos))
sc = {}
out = subprocess.run([TABIX, "-R", bed, SRC], capture_output=True, text=True).stdout
for l in out.splitlines():
    p = l.split("\t"); d = sc.setdefault((p[0], int(p[1])), {"ref": p[2], "alts": {}}); d["alts"][p[3]] = float(p[4])
cnt = collections.Counter()
with gzip.open(P + "/gpn_scores.v3.tmp.gz", "wt") as o:
    o.write(hdr)
    for r in rows:
        if r[6] == "ref_mismatch":
            kp = r[0].split(":"); kr, ka = kp[2], kp[3]; d = sc.get((r[1], int(r[2])))
            if len(kr) != 1 or len(ka) != 1: r[6] = "indel_no_snv_score"
            elif d and ka == d["ref"] and kr in d["alts"]: r[3] = f"{-d['alts'][kr]:.3f}"; r[6] = "swap_neg"
            elif d and rc(ka) == d["ref"] and rc(kr) in d["alts"]: r[3] = f"{-d['alts'][rc(kr)]:.3f}"; r[6] = "swap_rc_neg"
        cnt[r[6]] += 1
        o.write("\t".join(r) + "\n")
os.replace(P + "/gpn_scores.tsv.gz", P + "/gpn_scores.v2.tsv.gz")
os.replace(P + "/gpn_scores.v3.tmp.gz", P + "/gpn_scores.tsv.gz")
meta = json.load(open(P + "/gpn_scores.meta.json")); meta["v3_match_counts"] = dict(cnt); meta["v3_fix"] = "swap -> -score(alt=key_ref); indel -> NA"
json.dump(meta, open(P + "/gpn_scores.meta.json", "w"), indent=1)
open(P + "/gpn.done", "w").write(json.dumps(meta) + "\n")
print(json.dumps(dict(cnt)))
