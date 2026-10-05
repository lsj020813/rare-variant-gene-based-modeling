#!/usr/bin/env python
import os as _os
if not _os.environ.get("PROJECT_ROOT", "").strip():
    raise ValueError("PROJECT_ROOT must be set and nonblank")
import sys
N=sys.argv[1]; W=(_os.environ["PROJECT_ROOT"] + '/work/prs')
info={}
for l in open(W+"/ref/ldblk_1kg_eas/snpinfo_1kg_hm3"):
    p=l.split()
    if p[0]==N: info[int(p[2])]=(p[1],{p[3],p[4]})
o=open(f"{W}/geno_hm3/chr{N}.rs.bim","w"); m=open(f"{W}/geno_hm3/chr{N}.rsmap","w"); k=0
for l in open(f"{W}/geno_hm3/chr{N}.pvar"):
    if l.startswith("#"): continue
    c,pos,vid,ref,alt=l.split()[:5]; pos=int(pos)
    if pos in info and {ref,alt}==info[pos][1]:
        rs=info[pos][0]; o.write(f"{N}\t{rs}\t0\t{pos}\t{alt}\t{ref}\n"); m.write(f"{rs}\t{vid}\n"); k+=1
o.close(); m.close(); print(N,k)
