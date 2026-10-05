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


import gzip, json, os, glob, random, sys
import numpy as np
ROOT=_config_path("${PROJECT_ROOT}/work"); OUT=f"{ROOT}/fset/out"; G1=f"{ROOT}/gate1/out"
SEED=20260922; NTOT=1000; CAP=400; MINV=50
rng=random.Random(SEED)

ann={}
for N in range(1,23):
    with gzip.open(f"{OUT}/ann_chr{N}.tsv.gz","rt") as f:
        f.readline()
        for line in f:
            p=line.rstrip("\n").split("\t")
            if p[4]=="1": ann[p[0]]=(int(p[1]), int(p[2]), int(p[3]))
info={}
for N in range(1,23):
    with gzip.open(f"{OUT}/uni_chr{N}.tsv.gz","rt") as f:
        f.readline()
        for line in f:
            p=line.rstrip("\n").split("\t")
            if p[5]=="1": continue
            info[p[0]]=(N,int(p[1]),float(p[2]),float(p[3]),int(p[4]))
chrlen={}
domU={}; domA={}
for N in range(1,23):
    with gzip.open(f"{G1}/windows_chr{N}.json.gz","rt") as f: Wd=json.load(f)
    chrlen[N]=len(Wd)
    for g,ent in Wd.items():
        ks=sorted({rec[0] for b,lst in ent["bins"].items() for rec in lst if rec[0] in info})
        if not ks: continue
        kb=[k for k in ks if k in ann]
        if len(kb)>=MINV:
            domU[g]=ks; domA[g]=kb
forced=set()
gs=f"{G1}/gene_sample.csv"
if os.path.exists(gs):
    import csv
    with open(gs) as f:
        for row in csv.DictReader(f):
            g=row.get("gene")
            if g in domA: forced.add(g)
genes=sorted(domA)
nb=np.array([len(domA[g]) for g in genes])
chrs=np.array([info[domA[g][0]][0] for g in genes])
sizes=np.array([chrlen[cc] for cc in chrs])
ct=np.quantile(sizes,[1/3,2/3]); tert=np.digitize(sizes,ct)
nq=np.quantile(nb,[.25,.5,.75]); quart=np.digitize(nb,nq)
strata={}
for i,g in enumerate(genes): strata.setdefault((int(tert[i]),int(quart[i])),[]).append(g)
per=NTOT//len(strata)
pick=set(forced)
order=sorted(strata, key=lambda k:-len(strata[k]))
for k in order:
    pool=[g for g in strata[k] if g not in pick]
    have=sum(1 for g in strata[k] if g in pick)
    take=max(0, per-have)
    rng.shuffle(pool)
    pick.update(pool[:take])
while len(pick)<NTOT:
    added=False
    for k in order:
        pool=[g for g in strata[k] if g not in pick]
        if pool:
            rng.shuffle(pool); pick.add(pool[0]); added=True
            if len(pick)>=NTOT: break
    if not added: break
pickl=sorted(pick)
rows=[]
capped=0
for g in pickl:
    ks_all=domU[g]; ks_b=domA[g]
    sel_b=sorted(rng.sample(ks_b, CAP)) if len(ks_b)>CAP else ks_b
    if len(ks_b)>CAP: capped+=1
    extra=[k for k in ks_all if k not in set(sel_b)]
    n_extra=min(CAP, max(0, CAP-len(sel_b))+CAP)
    sel_a_extra=sorted(rng.sample(extra, min(len(extra), CAP))) if extra else []
    for k in sel_b:
        N,p38,maf,r2,ty=info[k]; a=ann[k]
        rows.append((k,N,p38,g,maf,r2,ty,1,1 if g in forced else 0,len(ks_b),len(ks_all),a[0],a[1],a[2]))
    for k in sel_a_extra:
        N,p38,maf,r2,ty=info[k]; a=ann.get(k,(0,0,0))
        rows.append((k,N,p38,g,maf,r2,ty,0,1 if g in forced else 0,len(ks_b),len(ks_all),a[0],a[1],a[2]))
hdr=["key","chr","pos38","domain","maf","r2","typed","in_ub","forced","domain_n_ub","domain_n_all",
     "has_ccre","n_tf","n_re2g"]
os.makedirs(f"{ROOT}/fset/primary", exist_ok=True)
with open(f"{ROOT}/fset/primary/primary_sample.tsv","w") as f:
    f.write("\t".join(hdr)+"\n")
    for r in rows: f.write("\t".join(map(str,r))+"\n")
meta={"seed":SEED,"n_strata":len(strata),"per_stratum":per,"eligible_domains":len(genes),
      "sampled_domains":len(pickl),"forced_included":len(forced),"capped_domains":capped,
      "variants_ub":sum(1 for r in rows if r[7]==1),"variants_extra_for_ua":sum(1 for r in rows if r[7]==0),
      "cap":CAP,"min_variants_ub":MINV,
      "median_domain_n_ub":float(np.median([len(domA[g]) for g in pickl])),
      "median_domain_n_all":float(np.median([len(domU[g]) for g in pickl]))}
json.dump(meta, open(f"{ROOT}/fset/primary/primary_sample.meta.json","w"), indent=1)
print(json.dumps(meta))
