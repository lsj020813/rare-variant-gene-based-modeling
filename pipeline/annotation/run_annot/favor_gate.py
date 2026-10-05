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


import gzip, json, os, random, subprocess, sys, time
R = _config_path("${PROJECT_ROOT}/work/ref"); T0 = time.time()
def log(*a): print(f"[{time.time()-T0:6.0f}s]", *a, flush=True)

obs = {}
for line in open(f"{R}/annot/map38/chr22.map.tsv"):
    k38 = line.rstrip("\n").split("\t")[1]
    ch, pos, ref, alt = k38.split(":")
    if len(ref) != 1 or len(alt) != 1: continue
    d = obs.setdefault(int(pos), {"ref": ref, "alts": set()})
    d["alts"].add(alt)
log(f"our chr22 SNV positions {len(obs):,}")

random.seed(20260906)
sample = sorted(random.sample(list(obs), 1000))
cand_unobs, cand_obs = {}, {}
for p in sample:
    ref, alts = obs[p]["ref"], obs[p]["alts"]
    for b in "ACGT":
        if b == ref: continue
        (cand_obs if b in alts else cand_unobs)[(p, ref, b)] = None
log(f"sample positions 1,000 | 관측 대립 후보 {len(cand_obs):,} | 미관측 대립 후보 {len(cand_unobs):,}")
posset = set(sample)

COLS = ["cons", "epi_active", "epi_repr", "epi_trans", "tf", "cage_prom", "genehancer", "linsight"]
IDX  = {"cons": 6, "epi_active": 9, "epi_repr": 10, "epi_trans": 11, "tf": 22, "cage_prom": 23,
        "genehancer": 32, "linsight": 33}
EXPECT = {1: "variant_vcf", 6: "apc_conservation", 9: "apc_epigenetics_active",
          10: "apc_epigenetics_repressed", 11: "apc_epigenetics_transcription",
          22: "apc_transcription_factor", 23: "cage_promoter", 32: "genehancer", 33: "linsight"}
proc = subprocess.Popen(["tar", "xzOf", f"{R}/favor/chr22.tar.gz"], stdout=subprocess.PIPE, bufsize=1 << 22)
hdr = proc.stdout.readline().decode().rstrip("\n").split(",")
log(f"header fields {len(hdr)}")
bad = {i: (hdr[i-1] if i <= len(hdr) else "MISSING") for i, name in EXPECT.items() if i > len(hdr) or hdr[i-1] != name}
if bad:
    print("GATE FAIL: header mismatch", json.dumps(bad, ensure_ascii=False)); sys.exit(3)
log("header verified")

found_unobs, found_obs = {}, {}
n = 0
for raw in proc.stdout:
    n += 1
    f = raw.decode("utf-8", "ignore").rstrip("\n").split(",")
    v = f[0]
    if not v.startswith("22-"): continue
    parts = v.split("-")
    if len(parts) != 4: continue
    _, p, ref, alt = parts
    if not p.isdigit() or int(p) not in posset: continue
    if len(ref) != 1 or len(alt) != 1: continue
    off = len(f) - 35
    vals = {}
    for cn, i in IDX.items():
        j = i + off if i >= 30 else i
        vals[cn] = f[j-1] if 0 < j <= len(f) else ""
    key = (int(p), ref, alt)
    if key in cand_unobs: found_unobs[key] = vals
    elif key in cand_obs: found_obs[key] = vals
    if n % 2_000_000 == 0: log(f"scanned {n:,} rows | hits unobs {len(found_unobs):,} obs {len(found_obs):,}")
proc.stdout.close(); proc.wait()
log(f"scan done rows {n:,}")

def nonnull(found):
    out = {}
    for cn in COLS:
        k = sum(1 for v in found.values() if v[cn].strip() not in ("", "NA", "."))
        out[cn] = dict(nonnull=k, pct=round(100 * k / max(len(found), 1), 1))
    return out

res = dict(
    positions=1000, rows_scanned=n,
    unobserved=dict(candidates=len(cand_unobs), found=len(found_unobs),
                    lookup_pct=round(100 * len(found_unobs) / max(len(cand_unobs), 1), 1), cols=nonnull(found_unobs)),
    observed=dict(candidates=len(cand_obs), found=len(found_obs),
                  lookup_pct=round(100 * len(found_obs) / max(len(cand_obs), 1), 1), cols=nonnull(found_obs)),
    sample_unobs=[{"variant": f"22-{p}-{r}-{a}", **{k: v[k] for k in ("cons", "gpn_placeholder") if k in v}}
                  for (p, r, a), v in list(found_unobs.items())[:3]],
)
json.dump(res, open(f"{R}/annot/favor_gate_chr22.json", "w"), indent=1)
log("RESULT", json.dumps({k: (v if not isinstance(v, dict) else {kk: vv for kk, vv in v.items() if kk != "cols"})
                          for k, v in res.items() if k != "sample_unobs"}, ensure_ascii=False))
log("unobs cols", json.dumps(res["unobserved"]["cols"], ensure_ascii=False))
log("obs   cols", json.dumps(res["observed"]["cols"], ensure_ascii=False))
print("FAVOR_GATE_DONE")
