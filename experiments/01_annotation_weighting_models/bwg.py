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


import gzip, os, sys, json, subprocess, bisect
from collections import defaultdict

CH = sys.argv[1]
OUT = _config_path("${PROJECT_ROOT}/work/ref15/groupfiles_bwg")
GTF = _config_path("${PROJECT_ROOT}/work/ref15/deductive/gencode.sorted.gtf.gz")
KEYED = _config_path(f'${{PROJECT_ROOT}}/work/ref15/lift38_keyed/chr{CH}.keyed38.vcf.gz')
BAND  = _config_path(f'${{PROJECT_ROOT}}/work/ref15/band_vcf/chr{CH}.band.vcf.gz')
RE2G = _config_path("${PROJECT_ROOT}/work/ref15/re2g")
BCF = "bcftools"
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
assert genes, f"GATE FAIL chr{CH}: no genes"
base2full = {g.split(".")[0]: g for g in genes}

band37 = set()
_bt = subprocess.run(f"{BCF} query -f '%CHROM:%POS:%REF:%ALT\\n' {BAND}",
                     shell=True, capture_output=True, text=True).stdout
for _l in _bt.split():
    band37.add(_l[3:] if _l.startswith("chr") else _l)
assert band37, f"GATE FAIL chr{CH}: band VCF empty"

txt = subprocess.run(f"{BCF} query -f '%ID\\t%CHROM\\t%POS\\n' {KEYED}",
                     shell=True, capture_output=True, text=True).stdout
var = []
for line in txt.splitlines():
    f = line.split("\t")
    if len(f) != 3 or f[1].replace("chr", "") != CH: continue
    k37 = f[0][3:] if f[0].startswith("chr") else f[0]
    if k37 not in band37: continue
    var.append((int(f[2]), k37))
var.sort()
assert var, f"GATE FAIL chr{CH}: no band variants after keyed join"
vpos = [p for p, _ in var]

links = defaultdict(list); nlink = 0
for fp in sorted(os.listdir(RE2G)):
    if not fp.endswith(".bed.gz"): continue
    with gzip.open(f"{RE2G}/{fp}", "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 14 or f[0] != f"chr{CH}": continue
            g = base2full.get(f[13].split(".")[0])
            if not g: continue
            links[g].append((int(f[1]), int(f[2]))); nlink += 1

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

amap = window_assign(3_000, True)
kept = {g: ks for g, ks in amap.items() if len(ks) >= 2}
path = f"{OUT}/chr{CH}.B_3kb_re2g.txt"
with open(path, "w") as fh:
    for g in sorted(kept):
        ks = sorted(kept[g])
        fh.write(f"{g} var " + " ".join(ks) + "\n")
        fh.write(f"{g} anno " + " ".join(["all"] * len(ks)) + "\n")

tot = sum(len(v) for v in kept.values())
uniq = len({k for v in kept.values() for k in v})
outside = {k for v in kept.values() for k in v} - band37
assert not outside, f"GATE FAIL chr{CH}: {len(outside)} variants outside band pool"
assert kept and uniq > 0, f"GATE FAIL chr{CH}: empty output"
summ = {"chrom": CH, "genes": len(kept), "pairs": tot, "distinct": uniq,
        "mean_m": round(tot/len(kept),1), "redundancy": round(tot/uniq,2)}
json.dump(summ, open(f"{OUT}/chr{CH}.summary.json","w"))
print(f"[chr{CH}] genes {len(kept):,} distinct {uniq:,} redundancy {summ['redundancy']}x")
print(f"B_WG_CHR{CH}_DONE")
