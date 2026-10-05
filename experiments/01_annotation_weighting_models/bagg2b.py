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
import sys, subprocess, collections, os, time
import numpy as np

HEART = os.environ.get("WD_HEARTBEAT", "")
DEADMAN_MAX_AGE = float(os.environ.get("WD_MAX_AGE", "60"))

CH = sys.argv[1]
R = _config_path("${PROJECT_ROOT}/work/ref")
BCF = "bcftools"
OUT = os.environ.get("BAGG_OUT", f"{R}/l1/B")
os.makedirs(OUT, exist_ok=True)
if os.path.exists(f"{OUT}/chr{CH}.B.done"):
    print(f"[chr{CH}] already done"); sys.exit(0)

def hwm_gb():
    for l in open("/proc/self/status"):
        if l.startswith("VmHWM"):
            return int(l.split()[1]) / 1048576
    return -1.0

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

rows = []
key2rows = collections.defaultdict(list)
key2genes = collections.defaultdict(list)
for gi, (g, keys) in enumerate(genes):
    pat = collections.defaultdict(list)
    for k in keys:
        p0 = k2p.get(k)
        pat[card.get(p0, "none") if p0 else "none"].append(k)
    assert sum(len(v) for v in pat.values()) == len(keys), f"GATE FAIL slots {g}"
    for kk in sorted(pat):
        ri = len(rows)
        rows.append((g, kk, len(pat[kk])))
        for k in pat[kk]:
            key2rows[k].append(ri)
            key2genes[k].append(gi)

NR = len(rows)
NG = len(genes)
print(f"[chr{CH}] genes {NG:,}  distinct keys {len(key2rows):,}  rows {NR:,}", flush=True)
print(f"[chr{CH}] 예상 메모리: 행렬 {NR*NS*4/1e9:.2f}G + 검산 {NG*NS*4/1e9:.2f}G", flush=True)

prog_every = max(1000, len(key2rows) // 20)

B = np.zeros((NR, NS), dtype=np.float32)
W = np.zeros((NG, NS), dtype=np.float32)

q = subprocess.Popen(f"{BCF} query -f '%CHROM:%POS:%REF:%ALT[\t%DS]\n' {R}/band_vcf/chr{CH}.band.vcf.gz",
                     shell=True, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
seen = set()
got = 0
for l in q.stdout:
    i = l.index("\t")
    k = l[:i]
    rr = key2rows.get(k)
    if rr is None or k in seen:
        continue
    ds = np.fromstring(l[i + 1:], dtype=np.float32, sep="\t")
    if ds.shape[0] != NS:
        raise SystemExit(f"GATE FAIL DS width {ds.shape[0]} != {NS} at {k}")
    for ri in rr:
        B[ri] += ds
    for gi in key2genes[k]:
        W[gi] += ds
    seen.add(k)
    got += 1
    if got % prog_every == 0:
        print(f"[chr{CH}] streamed {got:,}/{len(key2rows):,} ({got/len(key2rows)*100:.0f}%)  hwm={hwm_gb():.1f}G", flush=True)
        if HEART and os.path.exists(HEART):
            age = time.time() - os.path.getmtime(HEART)
            if age > DEADMAN_MAX_AGE:
                raise SystemExit(f"DEADMAN: watchdog heartbeat stale ({age:.0f}s > {DEADMAN_MAX_AGE}s) — aborting")
        elif HEART:
            raise SystemExit(f"DEADMAN: watchdog heartbeat missing ({HEART}) — aborting")
q.wait()

miss = len(key2rows) - len(seen)
assert miss == 0, f"GATE FAIL DS fetch: {miss} missing of {len(key2rows)}"
print(f"[chr{CH}] DS loaded {got:,}  hwm={hwm_gb():.1f}G", flush=True)

gi_rows = collections.defaultdict(list)
for ri, (g, kk, n) in enumerate(rows):
    gi_rows[g].append(ri)
gname2i = {g: i for i, (g, _) in enumerate(genes)}
part_fail = 0
for g, ris in gi_rows.items():
    psum = B[ris].sum(axis=0)
    if not np.allclose(W[gname2i[g]], psum, atol=1e-3):
        part_fail += 1
assert part_fail == 0, f"GATE FAIL partition identity: {part_fail} genes"
del W

LOCK = os.environ.get("BAGG_SAVE_LOCK", f"{OUT}/.save.lock")
import fcntl
_lk = open(LOCK, "w")
print(f"[chr{CH}] 저장 잠금 대기 (hwm={hwm_gb():.2f}G)", flush=True)
fcntl.flock(_lk, fcntl.LOCK_EX)
print(f"[chr{CH}] 저장 잠금 획득", flush=True)

M = B.astype(np.float16)
del B
np.savez_compressed(f"{OUT}/chr{CH}.B.npz",
                    B=M, rows=np.array([f"{g}|{k}|{n}" for g, k, n in rows]))
fcntl.flock(_lk, fcntl.LOCK_UN); _lk.close()
open(f"{OUT}/chr{CH}.B.done", "w").write(f"ok rows={NR}\n")
peak = hwm_gb()
with open(f"{OUT}/chr{CH}.peak_rss.txt", "w") as fh:
    fh.write(f"peak_rss_gb={peak:.2f}\nrows={NR}\nkeys={len(key2rows)}\ngenes={NG}\n")
print(f"[chr{CH}] B rows {NR:,}  ({M.nbytes/1e9:.2f} GB fp16)  GATES OK  peak_rss={peak:.2f}G", flush=True)
print(f"B_CHR{CH}_DONE")
