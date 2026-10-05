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


import numpy as np, glob, csv, math
W=_config_path("${PROJECT_ROOT}/work")
z=np.load(f"{W}/ref/annot/cache/fm_all.npz", allow_pickle=True)
print("npz fields:", list(z.keys()))
print("cols:", [str(c) for c in z["cols"]])
m=z["chr"]==19
key37=set(map(str,z["key37"][m]))
k38=None
for f in z.keys():
    if "38" in f: print("has field", f); k38=set(map(str,z[f][m]))
gk={}; 
for line in open(f"{W}/run_l3b/out/G_phi.txt"):
    t=line.rstrip("\n").split()
    if len(t)>2 and t[1]=="var": gk[t[0]]=t[2:]
allv=[v for vs in gk.values() for v in vs]
print("G_phi genes", len(gk), "vars", len(allv), "unique", len(set(allv)))
sz=set()
for f in glob.glob(f"{W}/run_l3b/out/arm_none/part*.singleAssoc.txt"):
    for r in csv.DictReader(open(f), delimiter="\t"):
        sz.add(f'{r["CHR"]}:{r["POS"]}:{r["Allele1"]}:{r["Allele2"]}')
print("singleAssoc keys", len(sz))
S=set(allv)
print("G_phi∩singleAssoc", len(S&sz), " G_phi∩key37", len(S&key37), " G_phi∩key38", (len(S&k38) if k38 else None))
ex=next(iter(S)); print("key style: chr-prefixed?", ex.startswith("chr"), " nfields", len(ex.split(":")))
ex2=next(iter(sz)); print("singleAssoc style: chr-prefixed?", ex2.startswith("chr"))
pos=np.array([int(k.split(":")[1]) for k in S]); print("G_phi pos min/max", pos.min(), pos.max())
pos2=np.array([int(k.split(":")[1]) for k in sz]); print("sA pos min/max", pos2.min(), pos2.max())
T="ENSG00000129353.15 ENSG00000213892.12 ENSG00000186567.14 ENSG00000130202.10 ENSG00000130204.13 ENSG00000104856.15 ENSG00000069399.15 ENSG00000079805.19 ENSG00000142453.13 ENSG00000127616.22 ENSG00000129354.12".split()
print("target genes present:", sum(g in gk for g in T), "/", len(T), " target vars:", sum(len(gk[g]) for g in T if g in gk))
rng=np.random.default_rng(20260910)
others=sorted(g for g in gk if g not in T)
ctrl=list(rng.choice(others, 20, replace=False))
print("ctrl20 vars:", sum(len(gk[g]) for g in ctrl))
print("ctrl20:", " ".join(ctrl))
print("ref fasta contigs:", open(f"{W}/tmp/ref/GRCh38_no_alt.fa.fai").read().count("\n"))
print("CHECK_DONE")
