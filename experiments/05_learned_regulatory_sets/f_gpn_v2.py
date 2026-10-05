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


import gzip, os, sys, json, subprocess, collections, time
ROOT = _config_path("${PROJECT_ROOT}/work"); OUT = f"{ROOT}/fset/primary"
SRC = f"{ROOT}/ref/gpnmsa/scores.tsv.bgz"
TABIX = "tabix"
RC = {"A": "T", "C": "G", "G": "C", "T": "A"}
def rc(s): return "".join(RC.get(b, "N") for b in reversed(s))
t0 = time.time()
rows = []
with open(sys.argv[1]) as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        if p[0] == "key" or len(p) < 3: continue
        kp = p[0].split(":")
        ref, alt = (kp[2], kp[3]) if len(kp) >= 4 else ("", "")
        rows.append((p[1], int(p[2]), p[0], ref, alt))
bypos = collections.defaultdict(list)
for ch, po, k, r, a in rows: bypos[(ch, po)].append((k, r, a))
def chrsort(c):
    return (0, int(c)) if c.isdigit() else (1, c)
bed = f"{OUT}/gpn_positions.bed"
with open(bed, "w") as o:
    for ch, po in sorted(bypos, key=lambda x: (chrsort(x[0]), x[1])):
        o.write(f"{ch}\t{po-1}\t{po}\n")
n_pos = len(bypos)
scores = {}
n_lines = 0
proc = subprocess.Popen([TABIX, "-R", bed, SRC], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
for line in proc.stdout:
    n_lines += 1
    p = line.rstrip("\n").split("\t")
    if len(p) < 5: continue
    kk = (p[0].replace("chr", ""), int(p[1]))
    if kk not in bypos: continue
    d = scores.setdefault(kk, {"ref": p[2], "alts": {}})
    try: d["alts"][p[3]] = float(p[4])
    except ValueError: pass
rc_code = proc.wait()
cnt = collections.Counter()
tmp = f"{OUT}/gpn_scores.tsv.gz.tmp"
with gzip.open(tmp, "wt") as o:
    o.write("key\tchr\tpos38\tgpn_score\tgpn_min3\tgpn_mean3\tgpn_match\n")
    for ch, po, k, r, a in rows:
        d = scores.get((ch, po))
        sc = ""; mn = ""; me = ""; m = "no_record"
        if d and d["alts"]:
            vals = list(d["alts"].values()); mn = f"{min(vals):.3f}"; me = f"{sum(vals)/len(vals):.3f}"
            if r == d["ref"] and a in d["alts"]: sc = f"{d['alts'][a]:.3f}"; m = "alt"
            elif rc(r) == d["ref"] and rc(a) in d["alts"]: sc = f"{d['alts'][rc(a)]:.3f}"; m = "alt_rc"
            else: m = "ref_mismatch"
        cnt[m] += 1
        o.write(f"{k}\t{ch}\t{po}\t{sc}\t{mn}\t{me}\t{m}\n")
os.replace(tmp, f"{OUT}/gpn_scores.tsv.gz")
meta = {"tabix_rc": rc_code, "records_read": n_lines, "positions_requested": n_pos, "positions_found": len(scores),
        "variants": len(rows), "match_counts": dict(cnt), "seconds": round(time.time() - t0, 1), "method": "tabix -R (index present)"}
json.dump(meta, open(f"{OUT}/gpn_scores.meta.json", "w"), indent=1)
print(json.dumps(meta))
with open(f"{OUT}/gpn.done", "w") as o: o.write(json.dumps(meta) + "\n")
