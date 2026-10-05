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

import glob, gzip, zipfile, os, json, collections, time
D = _config_path("${PROJECT_ROOT}/work/ref/bbj"); OUT = _config_path("${PROJECT_ROOT}/work/ref/annot/bbj_l1/raw"); os.makedirs(OUT, exist_ok=True)
WANT = _config_path("${PROJECT_ROOT}/work/ref/bbj_fm/union/keys37.txt")
T0 = time.time()
want = set(l.strip() for l in open(WANT) if l.strip()); assert len(want) == _config_number("N_BBJ_VARIANTS", int, True)
info = {}; n = collections.Counter(); diff = collections.defaultdict(float); orient = {}
per_trait = {}
for z in sorted(glob.glob(f"{D}/hum0197.v3.BBJ.*.v1.zip")):
    tr = z.split(".BBJ.")[1].split(".v1")[0]
    zf = zipfile.ZipFile(z); names = [x for x in zf.namelist() if x.endswith("auto.txt.gz")]
    assert len(names) == 1, f"GATE FAIL {tr} {zf.namelist()}"
    rows = hits = fwd = rev = bad = 0
    with gzip.open(zf.open(names[0]), "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        if "INFO" in hdr: ii, ic, ip, i1, i0, iN = [hdr.index(x) for x in ("INFO", "CHR", "BP", "ALLELE1", "ALLELE0")] + [None]; schema = "BOLT"
        else: ii, ic, ip, i1, i0, iN = [hdr.index(x) for x in ("imputationInfo", "CHR", "POS", "Allele1", "Allele2", "N")]; schema = "SAIGE"
        nmax = 0
        for line in fh:
            f = line.rstrip("\n").split("\t"); rows += 1
            if iN is not None:
                try: nmax = max(nmax, int(float(f[iN])))
                except ValueError: pass
            k1 = f"{f[ic]}:{f[ip]}:{f[i1]}:{f[i0]}"
            if k1 in want: k, o = k1, 0
            else:
                k2 = f"{f[ic]}:{f[ip]}:{f[i0]}:{f[i1]}"
                if k2 in want: k, o = k2, 1
                else: continue
            try: v = float(f[ii])
            except ValueError: bad += 1; continue
            hits += 1; fwd += (o == 0); rev += (o == 1)
            if k in info: diff[k] = max(diff[k], abs(info[k] - v))
            else: info[k] = v; orient[k] = o
            n[k] += 1
    per_trait[tr] = dict(schema=schema, rows=rows, hits=hits, fwd=fwd, rev=rev, nonnumeric=bad, N_max=(nmax if iN is not None else None))
    print(f"[rsq {time.time()-T0:5.0f}s] {tr} ({schema}, Nmax {nmax}): 행 {rows:,} 히트 {hits:,} (정 {fwd:,} 역 {rev:,}) 누적 고유 {len(info):,}", flush=True)
with open(f"{OUT}/variant_rsq.tsv.tmp", "w") as o:
    o.write("key37_bbj\tinfo\tn_traits\tmax_absdiff\torientation\n")
    for k, v in info.items(): o.write(f"{k}\t{v}\t{n[k]}\t{diff.get(k,0.0)}\t{orient[k]}\n")
os.replace(f"{OUT}/variant_rsq.tsv.tmp", f"{OUT}/variant_rsq.tsv")
incons = sum(1 for k in info if diff.get(k, 0) > 1e-6)
out_of_range = sum(1 for v in info.values() if not (0 <= v <= 1))
st = dict(union=len(want), covered=len(info), coverage_pct=round(len(info)/len(want)*100, 3), inconsistent_gt1e6=incons, out_of_range=out_of_range, per_trait=per_trait)
json.dump(st, open(f"{OUT}/rsq.stats.json", "w"), indent=1)
open(f"{OUT}/rsq.done", "w").write(f"ok {len(info)}\n")
print(json.dumps({k: v for k, v in st.items() if k != "per_trait"})); print("RSQ_DONE")
