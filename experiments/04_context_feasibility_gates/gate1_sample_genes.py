#!/usr/bin/env python
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import gzip, json, os, numpy as np, collections
ROOT=_config_path("${PROJECT_ROOT}/work"); OUT=f"{ROOT}/gate1/out"
SEED=20260921; PER=25; TARGET=300
rng=np.random.default_rng(SEED)
rows=[]
for N in range(1,23):
    wj=f"{OUT}/windows_chr{N}.json.gz"; vp=f"{OUT}/vpos_chr{N}.tsv.gz"
    if not (os.path.exists(wj) and os.path.exists(vp)): continue
    coding={}
    with gzip.open(vp,"rt") as f:
        f.readline()
        for line in f:
            a,b,cd=line.rstrip("\n").split("\t"); coding[a]=int(cd)
    with gzip.open(wj,"rt") as f: Wd=json.load(f)
    for g,ent in Wd.items():
        nb={}
        for b in ("B1","B2","B3","B4"):
            nb[b]=sum(1 for r in ent["bins"].get(b,[]) if coding.get(r[0],1)==0)
        if nb["B2"]>=2:
            rows.append({"gene":g,"chr":N,"nB1":nb["B1"],"nB2":nb["B2"],"nB3":nb["B3"],"nB4":nb["B4"],
                         "nU3":sum(nb.values())})
if not rows: raise SystemExit("no evaluable genes yet")
import pandas as pd
df=pd.DataFrame(rows)
chr_n=df.groupby("chr").size().sort_values(ascending=False)
tert={}
for i,ch in enumerate(chr_n.index): tert[ch]=i%3 if False else (0 if i<len(chr_n)/3 else (1 if i<2*len(chr_n)/3 else 2))
df["chr_tert"]=df["chr"].map(tert)
df["nB2_q"]=pd.qcut(df.nB2, 4, labels=False, duplicates="drop")
picks=[]
groups={k:v for k,v in df.groupby(["chr_tert","nB2_q"])}
deficit=0
for k in sorted(groups):
    sub=groups[k]
    take=min(PER, len(sub))
    idx=rng.choice(len(sub), size=take, replace=False)
    picks.append(sub.iloc[sorted(idx)])
    deficit += PER-take
if deficit>0:
    big=max(groups, key=lambda k: len(groups[k]))
    already=set(pd.concat(picks).gene)
    pool=groups[big][~groups[big].gene.isin(already)]
    take=min(deficit,len(pool))
    if take>0:
        idx=rng.choice(len(pool), size=take, replace=False)
        picks.append(pool.iloc[sorted(idx)])
S=pd.concat(picks).drop_duplicates("gene")
S.to_csv(f"{OUT}/gene_sample.csv", index=False)
for N,grp in S.groupby("chr"):
    with open(f"{OUT}/genes_chr{N}.txt","w") as f: f.write("\n".join(grp.gene)+"\n")
df.to_csv(f"{OUT}/evaluable_genes_all.csv", index=False)
print(json.dumps({"evaluable_genes":len(df),"sampled":len(S),"target":TARGET,"seed":SEED,
                  "per_stratum":PER,"strata":len(groups),
                  "chr_covered":int(S.chr.nunique()),
                  "median_nB2_sampled":float(S.nB2.median()),"median_nU3_sampled":float(S.nU3.median())}))
