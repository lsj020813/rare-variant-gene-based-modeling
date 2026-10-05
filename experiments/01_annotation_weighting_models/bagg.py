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
import sys, subprocess, collections, os
import numpy as np
CH = sys.argv[1]
R = _config_path("${PROJECT_ROOT}/work/ref")
BCF = "bcftools"
OUT = f"{R}/l1/B"
os.makedirs(OUT, exist_ok=True)
if os.path.exists(f"{OUT}/chr{CH}.B.done"):
    print(f"[chr{CH}] already done"); sys.exit(0)

card = {}
with open(f"{R}/l1/chr{CH}.ccre_card.tsv") as fh:
    fh.readline()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        card[f[1]] = f[2]
out = subprocess.run(f"{BCF} query -f '%POS\t%ID\n' {R}/lift38_keyed/chr{CH}.keyed38.vcf.gz",
                     shell=True, capture_output=True, text=True)
k2p = {}
for l in out.stdout.splitlines():
    p, kid = l.split("\t")
    k2p[kid[3:] if kid.startswith("chr") else kid] = str(int(p) - 1)

genes = []
for line in open(f"{R}/groupfiles_bwg/chr{CH}.B_3kb_re2g.txt"):
    f = line.split()
    if f[1] == "var":
        genes.append((f[0], [k[3:] if k.startswith("chr") else k for k in f[2:]]))

samples = subprocess.run(f"{BCF} query -l {R}/band_vcf/chr{CH}.band.vcf.gz",
                         shell=True, capture_output=True, text=True).stdout.split()
NS = len(samples)
assert NS == N_SAMPLES, f"GATE FAIL samples {NS}"

need = {}
for g, keys in genes:
    for k in keys: need.setdefault(k, None)
print(f"[chr{CH}] genes {len(genes):,}  distinct keys {len(need):,}", flush=True)
q = subprocess.Popen(f"{BCF} query -f '%CHROM:%POS:%REF:%ALT[\t%DS]\n' {R}/band_vcf/chr{CH}.band.vcf.gz",
                     shell=True, stdout=subprocess.PIPE, text=True, bufsize=1<<20)
got = 0
for l in q.stdout:
    i = l.index("\t")
    k = l[:i]
    if k in need and need[k] is None:
        need[k] = np.fromstring(l[i+1:], dtype=np.float32, sep="\t")
        got += 1
q.wait()
miss = sum(1 for v in need.values() if v is None)
assert miss == 0, f"GATE FAIL DS fetch: {miss} missing"
print(f"[chr{CH}] DS loaded {got:,}", flush=True)

rows = []; mats = []; part_fail = 0
for g, keys in genes:
    assert len(keys) > 0
    pat = collections.defaultdict(list)
    for k in keys:
        p0 = k2p.get(k)
        pat[card.get(p0, "none") if p0 else "none"].append(k)
    assert sum(len(v) for v in pat.values()) == len(keys), f"GATE FAIL slots {g}"
    whole = np.zeros(NS, dtype=np.float32)
    psum = np.zeros(NS, dtype=np.float32)
    for kk in sorted(pat):
        b = np.zeros(NS, dtype=np.float32)
        for k in pat[kk]: b += need[k]
        rows.append((g, kk, len(pat[kk])))
        mats.append(b.astype(np.float16))
        psum += b
    for k in keys: whole += need[k]
    if not np.allclose(whole, psum, atol=1e-3): part_fail += 1
assert part_fail == 0, f"GATE FAIL partition identity: {part_fail} genes"
M = np.stack(mats)
np.savez_compressed(f"{OUT}/chr{CH}.B.npz",
                    B=M, rows=np.array([f"{g}|{k}|{n}" for g, k, n in rows]))
open(f"{OUT}/chr{CH}.B.done", "w").write(f"ok rows={len(rows)}\n")
print(f"[chr{CH}] B rows {len(rows):,}  ({M.nbytes/1e9:.2f} GB fp16)  GATES OK", flush=True)
print(f"B_CHR{CH}_DONE")
