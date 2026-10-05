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


import subprocess, json, random, time, collections
R = _config_path("${PROJECT_ROOT}/work/ref"); T0 = time.time()
def log(*a): print(f"[{time.time()-T0:6.0f}s]", *a, flush=True)
obs = {}
for line in open(f"{R}/annot/map38/chr22.map.tsv"):
    k38 = line.rstrip("\n").split("\t")[1]
    ch, pos, ref, alt = k38.split(":")
    if len(ref) == 1 and len(alt) == 1:
        obs.setdefault(int(pos), {"ref": ref, "alts": set()})["alts"].add(alt)
random.seed(20260906)
sample = sorted(random.sample(list(obs), 1000)); posset = set(sample)
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
proc = subprocess.Popen(["tar", "xzOf", f"{R}/favor/chr22.tar.gz"], stdout=subprocess.PIPE, bufsize=1 << 22)
proc.stdout.readline()
seen = collections.defaultdict(set)
n = 0
for raw in proc.stdout:
    n += 1
    v = raw[:40].decode("utf-8", "ignore").split(",", 1)[0]
    if not v.startswith("22-"): continue
    p = v.split("-")
    if len(p) != 4 or not p[1].isdigit(): continue
    if int(p[1]) in posset and len(p[2]) == 1 and len(p[3]) == 1:
        seen[int(p[1])].add((p[2], p[3]))
proc.stdout.close(); proc.wait(); log(f"scan done rows {n:,}")
pos_found = len(seen)
ref_exact = ref_comp = ref_other = 0
alts_per_pos = collections.Counter()
for p in seen:
    refs = {r for r, a in seen[p]}
    alts_per_pos[len(seen[p])] += 1
    our = obs[p]["ref"]
    if our in refs: ref_exact += 1
    elif COMP[our] in refs: ref_comp += 1
    else: ref_other += 1
res = dict(rows=n, sampled=1000, pos_found=pos_found, pos_found_pct=round(100*pos_found/1000, 1),
           ref_exact=ref_exact, ref_complement=ref_comp, ref_other=ref_other,
           alleles_per_found_pos=dict(sorted(alts_per_pos.items())),
           verdict=("H2 대립불일치" if pos_found > 900 else "H1 커버리지" if pos_found < 400 else "혼합"))
json.dump(res, open(f"{R}/annot/favor_gate2_chr22.json", "w"), indent=1)
log("RESULT", json.dumps(res, ensure_ascii=False))
print("GATE2_DONE")
