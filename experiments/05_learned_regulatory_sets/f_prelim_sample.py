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
ROOT=_config_path("${PROJECT_ROOT}/work"); OUT=f"{ROOT}/fset/prelim"; G1=f"{ROOT}/gate1/out"; UNI=f"{ROOT}/fset/out"
os.makedirs(OUT, exist_ok=True)
SEED=20260922; NDOM=150; CAP=200; MINV=50
rng=random.Random(SEED)
chroms=sorted(int(f.split("uni_chr")[1].split(".tsv.gz")[0]) for f in glob.glob(f"{UNI}/uni_chr*.tsv.gz"))
info={}
for N in chroms:
    with gzip.open(f"{UNI}/uni_chr{N}.tsv.gz","rt") as f:
        f.readline()
        for line in f:
            p=line.rstrip("\n").split("\t")
            if p[5]=="1": continue
            info[p[0]]=(N,int(p[1]),float(p[2]),float(p[3]),int(p[4]))
dom={}
for N in chroms:
    with gzip.open(f"{G1}/windows_chr{N}.json.gz","rt") as f: Wd=json.load(f)
    for g,ent in Wd.items():
        ks=[rec[0] for b,lst in ent["bins"].items() for rec in lst if rec[0] in info]
        if len(ks)>=MINV: dom[g]=sorted(set(ks))
pool=sorted(dom)
pick=rng.sample(pool, min(NDOM,len(pool)))
rows=[]
for g in pick:
    ks=dom[g]
    sel=sorted(rng.sample(ks, CAP)) if len(ks)>CAP else ks
    for k in sel:
        N,p38,maf,r2,ty=info[k]
        rows.append((k,N,p38,g,maf,r2,ty,len(ks)))
with open(f"{OUT}/prelim_sample.tsv","w") as f:
    f.write("key\tchr\tpos38\tdomain\tmaf\tr2\ttyped\tdomain_n_total\n")
    for r in rows: f.write("\t".join(map(str,r))+"\n")
json.dump({"seed":SEED,"chroms":chroms,"eligible_domains":len(pool),"sampled_domains":len(pick),
           "variants":len(rows),"cap":CAP,"min_variants":MINV,
           "capped_domains":sum(1 for g in pick if len(dom[g])>CAP)},
          open(f"{OUT}/prelim_sample.meta.json","w"))
print(json.dumps(json.load(open(f"{OUT}/prelim_sample.meta.json"))))
