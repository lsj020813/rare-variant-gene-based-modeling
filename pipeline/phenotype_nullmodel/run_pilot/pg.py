#!/usr/bin/env python3
import os
import json

def required(name):
    value = os.environ.get(name, "")
    if not value.strip():
        raise RuntimeError("Required environment variable is missing or empty: " + name)
    return value

PROJECT_ROOT = required("PROJECT_ROOT")
import gzip, os, sys, json, subprocess, bisect
from collections import defaultdict

CH = "19"
OUT = f"{PROJECT_ROOT}/work/ref/groupfiles_pilot"
GTF = f"{PROJECT_ROOT}/work/ref/deductive/gencode.sorted.gtf.gz"
KEYED = f"{PROJECT_ROOT}/work/ref/lift38_keyed/chr{CH}.keyed38.vcf.gz"
BAND  = f"{PROJECT_ROOT}/work/ref/band_vcf/chr{CH}.band.vcf.gz"
RE2G = f"{PROJECT_ROOT}/work/ref/re2g"
BCF = os.environ.get("BCFTOOLS", "bcftools")
SMOKE = "--smoke" in sys.argv
os.makedirs(OUT, exist_ok=True)

genes = {}
with gzip.open(GTF, "rt") as fh:
    for line in fh:
        if line[0] == "#": continue
        f = line.split("\t")
        if f[2] != "gene" or f[0] != f"chr{CH}": continue
        if 'gene_type "protein_coding"' not in f[8]: continue
        i = f[8].find('gene_id "'); gid = f[8][i+9:f[8].find('"', i+9)]
        s, e = int(f[3]), int(f[4])
        genes[gid] = {"tss": s if f[6] == "+" else e, "start": s, "end": e}
print(f"protein-coding genes chr{CH}: {len(genes):,}", flush=True)
assert genes, "GATE FAIL: no genes parsed"

base2full = {g.split(".")[0]: g for g in genes}
tss_sorted = sorted((v["tss"], g) for g, v in genes.items())
tss_pos = [t for t, _ in tss_sorted]

band37 = set()
_bt = subprocess.run([BCF, "query", "-f", "%CHROM:%POS:%REF:%ALT\n", BAND],
                     check=True, capture_output=True, text=True).stdout
for _l in _bt.split():
    band37.add(_l[3:] if _l.startswith("chr") else _l)
print(f"band VCF pool (what V1 reads): {len(band37):,}", flush=True)
assert band37, "GATE FAIL: band VCF read produced nothing"

txt = subprocess.run(
    [BCF, "query", "-f", "%ID\t%CHROM\t%POS\n", KEYED], check=True,
    capture_output=True, text=True).stdout
var = []
for line in txt.splitlines():
    f = line.split("\t")
    if len(f) != 3 or f[1].replace("chr", "") != CH: continue
    k37 = f[0][3:] if f[0].startswith("chr") else f[0]
    if k37 not in band37: continue
    var.append((int(f[2]), k37))
var.sort()
print(f"band variants chr{CH} (38-mapped): {len(var):,}", flush=True)
assert var, "GATE FAIL: no band variants"
if SMOKE: var = var[:20000]

vpos = [p for p, _ in var]

links = defaultdict(list)
nlink = 0
for fp in sorted(os.listdir(RE2G)):
    if not fp.endswith(".bed.gz"): continue
    with gzip.open(f"{RE2G}/{fp}", "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 14 or f[0] != f"chr{CH}": continue
            g = base2full.get(f[13].split(".")[0])
            if not g: continue
            links[g].append((int(f[1]), int(f[2]))); nlink += 1
print(f"rE2G links on chr{CH} mapped to pc genes: {nlink:,} over {len(links):,} genes", flush=True)
assert nlink > 0, "GATE FAIL: rE2G join produced nothing"

def window_assign(halfwidth, use_re2g):
    out = defaultdict(set)
    for g, gi in genes.items():
        lo = min(gi["tss"] - halfwidth, gi["start"])
        hi = max(gi["tss"] + halfwidth, gi["end"])
        i, j = bisect.bisect_left(vpos, lo), bisect.bisect_right(vpos, hi)
        for p, k in var[i:j]: out[g].add(k)
        if use_re2g:
            for es, ee in links.get(g, ()):
                a, b = bisect.bisect_left(vpos, es), bisect.bisect_right(vpos, ee)
                for p, k in var[a:b]: out[g].add(k)
    return out

def nearest_assign():
    out = defaultdict(set)
    for p, k in var:
        i = bisect.bisect_left(tss_pos, p)
        best, bd = None, None
        for j in (i-1, i):
            if 0 <= j < len(tss_sorted):
                d = abs(tss_sorted[j][0] - p)
                if bd is None or d < bd: bd, best = d, tss_sorted[j][1]
        if best: out[best].add(k)
    return out

ARMS = {
    "B_3kb_re2g": lambda: window_assign(3_000, True),
    "C_nearest":  nearest_assign,
    "D_3kb_only": lambda: window_assign(3_000, False),
}

summary = {}
for name, fn in ARMS.items():
    amap = fn()
    kept = {g: ks for g, ks in amap.items() if len(ks) >= 2}
    sfx = ".smoke" if SMOKE else ""
    path = f"{OUT}/chr{CH}.{name}{sfx}.txt"
    with open(path, "w") as fh:
        for g in sorted(kept):
            ks = sorted(kept[g])
            fh.write(f"{g} var " + " ".join(ks) + "\n")
            fh.write(f"{g} anno " + " ".join(["all"] * len(ks)) + "\n")
    tot = sum(len(v) for v in kept.values())
    uniq = len({k for v in kept.values() for k in v})
    summary[name] = {"genes": len(kept), "gene_variant_pairs": tot,
                     "distinct_variants": uniq,
                     "mean_m": round(tot/len(kept), 1) if kept else 0,
                     "redundancy": round(tot/uniq, 2) if uniq else 0}
    print(f"[{name}] genes {len(kept):,} | pairs {tot:,} | distinct {uniq:,} | "
          f"mean m {summary[name]['mean_m']} | redundancy {summary[name]['redundancy']}x", flush=True)

json.dump(summary, open(f"{OUT}/summary{'_smoke' if SMOKE else ''}.json", "w"), indent=1)

import glob as _g
A_keys = set()
for _fp in _g.glob(f"{PROJECT_ROOT}/work/ref/groupfiles_chunks/chr%s.part*.txt" % CH):
    for _line in open(_fp):
        _f = _line.rstrip("\n").split()
        if len(_f) >= 3 and _f[1] == "var": A_keys.update(_f[2:])
allpilot = set()
for _n in ARMS:
    for _line in open(f"{OUT}/chr{CH}.{_n}{'.smoke' if SMOKE else ''}.txt"):
        _f = _line.rstrip("\n").split()
        if len(_f) >= 3 and _f[1] == "var": allpilot.update(_f[2:])
outside = allpilot - band37
summary["_baseline"] = {"A_distinct": len(A_keys), "pilot_union": len(allpilot),
                        "band_pool": len(band37), "pilot_outside_pool": len(outside),
                        "A_outside_pool": len(A_keys - band37)}
print(f"[baseline] A distinct {len(A_keys):,} | pilot union {len(allpilot):,} | "
      f"band pool {len(band37):,} | pilot outside pool {len(outside):,} | "
      f"A outside pool {len(A_keys - band37):,}", flush=True)

fail = []
if outside:
    fail.append(f"{len(outside):,} pilot variants outside the band pool — "
                "arms would answer a different exam than V1")
if not SMOKE and (A_keys - band37):
    fail.append(f"{len(A_keys - band37):,} V1 variants outside the band pool — "
                "baseline assumption broken; investigate before comparing")
for n in ARMS:
    a = summary[n]
    if a["genes"] == 0: fail.append(f"{n}: zero genes")
    if a["distinct_variants"] == 0: fail.append(f"{n}: zero variants — key join broken")
if summary["C_nearest"]["redundancy"] > 1.001:
    fail.append(f"C_nearest redundancy {summary['C_nearest']['redundancy']} > 1 — "
                "nearest assignment must be exclusive")
if summary["B_3kb_re2g"]["distinct_variants"] <= summary["D_3kb_only"]["distinct_variants"]:
    fail.append("B (3kb U rE2G) must cover MORE variants than D (3kb only)")
if fail:
    for m in fail: print("GATE FAIL:", m)
    sys.exit(6)
print("SMOKE_OK" if SMOKE else "PILOT_GROUPFILES_COMPLETE")
