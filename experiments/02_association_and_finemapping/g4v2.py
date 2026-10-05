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


import sys, os, json, subprocess, glob
from collections import defaultdict

SMOKE = "--smoke" in sys.argv
CH = ["21","22"] if SMOKE else [str(i) for i in range(1,23)]
THR = 0.1
IDX  = _config_path("${PROJECT_ROOT}/work/ref/gate4/pangolin_idx")
OUT  = _config_path("${PROJECT_ROOT}/work/ref/gate4")
GRP  = _config_path("${PROJECT_ROOT}/work/ref/groupfiles_chunks")
SAI  = (_config_path("${PROJECT_ROOT}/"
        "resources/annotation_data/spliceAI/spliceai_scores.raw.snv.hg38.vcf.gz"))
BCF  = "bcftools"

def pangolin(n):
    d = {}
    fp = f"{IDX}/chr{n}.tsv"
    if not os.path.exists(fp): return d
    for line in open(fp):
        f = line.rstrip("\n").split("\t")
        if len(f) == 3:
            d[f[0]] = (int(f[1]), float(f[2]))
    return d

def spliceai(n, keyset):
    d = {}
    p = subprocess.Popen(
        f"{BCF} query -r {n} -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/SpliceAI\\n' {SAI} 2>/dev/null",
        shell=True, stdout=subprocess.PIPE, text=True, bufsize=1<<20)
    for line in p.stdout:
        f = line.rstrip("\n").split("\t")
        if len(f) < 5: continue
        k = f"{f[0].replace('chr','')}:{f[1]}:{f[2]}:{f[3]}"
        if k not in keyset: continue
        bv, bd = 0.0, 0
        for ann in f[4].split(","):
            q = ann.split("|")
            if len(q) < 6: continue
            try: ag, al, dg, dl = (float(x) if x not in ("",".") else 0.0 for x in q[2:6])
            except ValueError: continue
            for v, dd in ((ag,1),(dg,1),(al,-1),(dl,-1)):
                if v > bv: bv, bd = v, dd
        if bv >= THR: d[k] = (bd, bv)
    p.wait()
    return d

def gene_map(n):
    m = defaultdict(set)
    for fp in glob.glob(f"{GRP}/chr{n}.part*.txt"):
        for line in open(fp):
            f = line.rstrip("\n").split()
            if len(f) < 3 or f[1] != "var": continue
            for v in f[2:]:
                m[v].add(f[0])
    return m

def key37to38(n):
    txt = subprocess.run(
        _config_path(f"{BCF} query -f '%ID\\t%CHROM\\t%POS\\t%REF\\t%ALT\\n' ${{PROJECT_ROOT}}/work/ref/lift38_keyed/chr{n}.keyed38.vcf.gz"),
        shell=True, capture_output=True, text=True).stdout
    d = {}
    for line in txt.splitlines():
        f = line.split("\t")
        if len(f) == 5:
            k37 = f[0][3:] if f[0].startswith("chr") else f[0]
            d[k37] = f"{f[1].replace('chr','')}:{f[2]}:{f[3]}:{f[4]}"
    return d

res, agree_genes = {}, defaultdict(int)
for n in CH:
    pg = pangolin(n)
    sa = spliceai(n, set(pg))
    both = set(pg) & set(sa)
    agree = {k for k in both if pg[k][0] == sa[k][0]}
    gm, k38 = gene_map(n), key37to38(n)
    v38_to_genes = defaultdict(set)
    for k37, genes in gm.items():
        kk = k37[3:] if k37.startswith("chr") else k37
        k = k38.get(kk)
        if k: v38_to_genes[k] |= genes
    if gm and not v38_to_genes:
        print(f"GATE FAIL chr{n}: group-file keys ({next(iter(gm))}) matched no keyed38 ID"); sys.exit(6)
    gcount = defaultdict(int)
    for k in agree:
        for g in v38_to_genes.get(k, ()): gcount[g] += 1
    ge2 = sum(1 for g,c in gcount.items() if c >= 2)
    for g,c in gcount.items(): agree_genes[g] += c
    res[f"chr{n}"] = {"pangolin": len(pg), "spliceai": len(sa), "both": len(both),
                      "agree": len(agree), "genes_ge1": len(gcount), "genes_ge2": ge2}
    print(f"[chr{n}] pang {len(pg):,} | sai {len(sa):,} | both {len(both):,} | "
          f"agree {len(agree):,} | genes>=1 {len(gcount):,} >=2 {ge2:,}", flush=True)

T = {k: sum(v[k] for v in res.values()) for k in ("pangolin","spliceai","both","agree")}
T["genes_ge2_global"] = sum(1 for g,c in agree_genes.items() if c >= 2)
T["genes_ge1_global"] = len(agree_genes)
res["TOTAL"] = T
sfx = "_smoke" if SMOKE else ""
json.dump(res, open(f"{OUT}/gate4{sfx}.json","w"), indent=1)

if T["pangolin"] == 0: print("GATE FAIL: Pangolin index empty"); sys.exit(6)
if T["spliceai"] == 0: print("GATE FAIL: SpliceAI produced zero directional calls — key mismatch"); sys.exit(6)
if T["both"] == 0:     print("GATE FAIL: zero overlap between the two sources"); sys.exit(6)
if T["genes_ge1_global"] == 0: print("GATE FAIL: gene mapping produced nothing — 37/38 key join broken"); sys.exit(6)

rate = T["agree"]/T["both"]
print(f"\nTOTAL both-scored {T['both']:,} | agree {T['agree']:,} | rate {rate:.4f}")
print(f"gene reach: >=1 {T['genes_ge1_global']:,} | >=2 {T['genes_ge2_global']:,}")
print("SMOKE_OK" if SMOKE else "GATE4_COMPLETE")
