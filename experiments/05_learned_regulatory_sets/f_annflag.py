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


import gzip, json, os, sys, bisect
CHR=sys.argv[1]; ROOT=_config_path("${PROJECT_ROOT}/work"); OUT=f"{ROOT}/fset/out"
chrname=f"chr{CHR}"
rows=[]
with gzip.open(f"{OUT}/uni_chr{CHR}.tsv.gz","rt") as f:
    f.readline()
    for line in f:
        p=line.rstrip("\n").split("\t")
        if p[5]=="1": continue
        rows.append((int(p[1]), p[0]))
rows.sort()
P=[r[0] for r in rows]

def load(path, ci, si, ei, gz=True, skip=None):
    iv=[]
    op=gzip.open(path,"rt") if gz else open(path)
    with op as f:
        for line in f:
            if skip and line.startswith(skip): continue
            q=line.rstrip("\n").split("\t")
            if len(q)<=max(ci,si,ei) or q[ci]!=chrname: continue
            iv.append((int(q[si]), int(q[ei])))
    iv.sort(); return iv

def sweep(iv, positions):
    out=[0]*len(positions)
    if not iv: return out
    import heapq
    active=[]; j=0
    for i,p in enumerate(positions):
        p0=p-1
        while j<len(iv) and iv[j][0]<=p0:
            heapq.heappush(active, iv[j][1]); j+=1
        while active and active[0]<=p0:
            heapq.heappop(active)
        out[i]=len(active)
    return out

ccre=sweep(load(f"{ROOT}/ref/b6_cards/ccre.s.bed",0,1,2,gz=False), P)
rm=f"{OUT}/remap_by_chr/{chrname}.bed"
tf=sweep(load(rm,0,1,2,gz=False), P) if os.path.exists(rm) else [0]*len(P)
re2g=[0]*len(P)
D=f"{ROOT}/ref/re2g"
for fn in sorted(os.listdir(D)):
    if not fn.endswith(".bed.gz"): continue
    s=sweep(load(os.path.join(D,fn),0,1,2), P)
    re2g=[a+(1 if b else 0) for a,b in zip(re2g,s)]
tmp=f"{OUT}/ann_chr{CHR}.tsv.gz.tmp"
nb=0
with gzip.open(tmp,"wt") as out:
    out.write("key\thas_ccre\tn_tf\tn_re2g_biosamples\tann_bearing\n")
    for (p38,key),a,b,d in zip(rows,ccre,tf,re2g):
        ab=1 if (a>0 or b>0 or d>0) else 0
        nb+=ab
        out.write(f"{key}\t{1 if a else 0}\t{b}\t{d}\t{ab}\n")
os.replace(tmp,f"{OUT}/ann_chr{CHR}.tsv.gz")
print(json.dumps({"chr":CHR,"noncoding_assigned":len(rows),"ann_bearing":nb,
                  "frac_ann_bearing":round(nb/max(len(rows),1),4),
                  "frac_ccre":round(sum(1 for x in ccre if x)/max(len(rows),1),4),
                  "frac_tf":round(sum(1 for x in tf if x)/max(len(rows),1),4),
                  "frac_re2g":round(sum(1 for x in re2g if x)/max(len(rows),1),4)}))
