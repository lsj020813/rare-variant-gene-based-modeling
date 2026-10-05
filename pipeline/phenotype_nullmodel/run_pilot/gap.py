import os
import json

def required(name):
    value = os.environ.get(name, "")
    if not value.strip():
        raise RuntimeError("Required environment variable is missing or empty: " + name)
    return value

PROJECT_ROOT = required("PROJECT_ROOT")
import glob
G=f'{PROJECT_ROOT}/work/ref/groupfiles_chunks'
per={}
for fp in glob.glob(f'{G}/chr19.part*.txt'):
    for line in open(fp):
        f=line.rstrip("\n").split()
        if len(f)>=3 and f[1]=="var": per[f[0]]=per.get(f[0],0)+len(f)-2
import statistics as st
v=sorted(per.values())
print("A genes:", len(per), "| median m:", v[len(v)//2], "| p90:", v[int(len(v)*0.9)], "| max:", v[-1])
allv=set()
for fp in glob.glob(f'{G}/chr19.part*.txt'):
    for line in open(fp):
        f=line.rstrip("\n").split()
        if len(f)>=3 and f[1]=="var": allv.update(f[2:])
pos=sorted(int(k.split(":")[1]) for k in allv)
print("A variant span:", f"{pos[0]:,}", "-", f"{pos[-1]:,}")
gaps=[(pos[i+1]-pos[i]) for i in range(len(pos)-1)]
big=[g for g in gaps if g>100000]
print("gaps >100kb:", len(big), "| total gap bp:", f"{sum(big):,}", "| chr19 length ~58.6Mb")
