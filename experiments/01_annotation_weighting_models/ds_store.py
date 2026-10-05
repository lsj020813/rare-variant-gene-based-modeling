#!/usr/bin/env python3
import os as _cfg_os
import re as _cfg_re
import math as _cfg_math

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)

def _config_number(name, cast=int, positive=False):
    value = cast(_cfg_os.environ[name])
    if not _cfg_math.isfinite(value) or (positive and value <= 0):
        raise ValueError(name + " must be finite" + (" and positive" if positive else ""))
    return value

import os as _os
N_SAMPLES = _config_number("N_SAMPLES", int, True)
import sys, os, subprocess, json, time, numpy as np
N = sys.argv[1]; R = _config_path("${PROJECT_ROOT}/work/ref15")
BCF = "bcftools"
V = f"{R}/band_vcf/chr{N}.band.vcf.gz"; FM = f"{R}/annot/fm/chr{N}.fm.tsv"
OUT = f"{R}/annot/ds"; os.makedirs(OUT, exist_ok=True); O = f"{OUT}/chr{N}.ds.npz"
THR = 0.05
def require(c, m):
    if not c: raise RuntimeError(f"GATE FAIL: {m}")
want = set()
with open(FM) as f:
    next(f)
    for line in f: want.add(line.split("\t", 1)[0])
require(want, "no assigned variants")
samples = subprocess.run(f"{BCF} query -l {V}", shell=True, capture_output=True, text=True).stdout.split()
require(len(samples) == N_SAMPLES, f"sample count {len(samples)}")
t0 = time.time()
p = subprocess.Popen(f"{BCF} query -f '%CHROM:%POS:%REF:%ALT\t[%DS\t]\n' {V}", shell=True, stdout=subprocess.PIPE, text=True, bufsize=1<<20)
data, idx, ptr, keys = [], [], [0], []
seen = 0; nz_total = 0
for line in p.stdout:
    k, rest = line.split("\t", 1)
    if k.startswith("chr"): k = k[3:]
    if k not in want: continue
    a = np.array(rest.rstrip("\t\n").split("\t"), dtype=np.float32)
    require(a.shape[0] == N_SAMPLES, f"row width {a.shape[0]} at {k}")
    require(np.nanmin(a) >= 0 and np.nanmax(a) <= 2.0, f"DS out of range at {k}")
    nzi = np.nonzero(a > THR)[0]
    data.append(a[nzi]); idx.append(nzi.astype(np.int32)); ptr.append(ptr[-1] + nzi.size); keys.append(k)
    seen += 1; nz_total += nzi.size
    if seen % 5000 == 0: print(f"[chr{N}] {seen}/{len(want)} nz={nz_total} {time.time()-t0:.0f}s", flush=True)
p.wait(); require(p.returncode == 0, "bcftools exit")
require(seen == len(want), f"rows {seen} != assigned {len(want)}")
np.savez(O + ".tmp.npz", data=np.concatenate(data), indices=np.concatenate(idx), indptr=np.asarray(ptr, dtype=np.int64),
         keys=np.array(keys), samples=np.array(samples), thr=np.float32(THR))
os.replace(O + ".tmp.npz", O)
st = {"chr": N, "variants": seen, "samples": N_SAMPLES, "nonzero": nz_total, "nz_pct": round(nz_total / (seen * N_SAMPLES) * 100, 3),
      "bytes": os.path.getsize(O), "sec": round(time.time() - t0)}
json.dump(st, open(O + ".stats.json", "w")); open(O + ".done", "w").write(f"ok {seen}\n")
print("DSSTATS", json.dumps(st)); print(f"[chr{N}] DS_DONE")
