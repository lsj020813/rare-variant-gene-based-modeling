#!/usr/bin/env python
import os as _os
if not _os.environ.get("PROJECT_ROOT", "").strip():
    raise ValueError("PROJECT_ROOT must be set and nonblank")
import sys, os, subprocess, gzip, json, numpy as np
L=sys.argv[1]; W=(_os.environ["PROJECT_ROOT"] + '/work/prs'); B=_os.environ.get("BCFTOOLS_BIN", "bcftools")
SEED=20260926; CAP=1024
loc={l.split()[3]:l.split()[:3] for l in open(W+"/out/loci.bed")}
ch,s,e=loc[L]; s=int(s); e=int(e)
sp=[l.rstrip("\n").split("\t") for l in open(W+"/private/split.tsv")][1:]
ids=[a for a,_ in sp]; istrain=np.array([b=="train" for _,b in sp])
os.makedirs(W+"/private/loci",exist_ok=True); os.makedirs(W+"/out/loci_meta",exist_ok=True)
sf=f"{W}/private/loci/{L}.samples"; open(sf,"w").write("\n".join(ids)+"\n")
cmd=[B,"query","-r",f"{ch}:{s}-{e}","-S",sf,"-i","INFO/R2>=0.3 && INFO/MAF>=0.0008","-f","%POS\t%REF\t%ALT\t%INFO/MAF\t%INFO/R2\t%INFO/TYPED[\t%DS]\n",f"{_os.environ['PROJECT_ROOT']}/work/ref/orig_index/chr{ch}.vcf.gz"]
p=subprocess.Popen(cmd,stdout=subprocess.PIPE,text=True)
meta=[]; cols=[]
for line in p.stdout:
    f=line.rstrip("\n").split("\t")
    meta.append((int(f[0]),f[1],f[2],float(f[3]),float(f[4]),f[5]))
    cols.append(np.array(f[6:],dtype=np.float32))
p.wait(); os.remove(sf)
n0=len(meta)
if n0==0: print(json.dumps({"locus":L,"n":0})); sys.exit(0)
D=np.vstack(cols).T
af=D[istrain].mean(0)/2; mafT=np.minimum(af,1-af)
keys=[f"{ch}:{m[0]}:{m[1]}:{m[2]}" for m in meta]
k2p={}
want=set(keys)
with gzip.open(f"{_os.environ['PROJECT_ROOT']}/work/ref/lift38_keyed/chr{ch}.keyed38.vcf.gz","rt") as g:
    for l in g:
        if l[0]=="#": continue
        a=l.split("\t",3)
        kk=a[2][3:] if a[2].startswith("chr") else a[2]
        if kk in want and a[0].replace("chr","")==ch: k2p[kk]=int(a[1])
cds=[]
with gzip.open((_os.environ["PROJECT_ROOT"] + '/work/ref/deductive/gencode.nochr.gtf.gz'),"rt") as g:
    for l in g:
        if l[0]=="#": continue
        a=l.split("\t",5)
        if a[0]==ch and a[2]=="CDS": cds.append((int(a[3]),int(a[4])))
cds.sort(); import bisect; st=[x[0] for x in cds]
def incds(p):
    i=bisect.bisect_right(st,p)-1
    while i>=0 and i>=len(st)-10**9:
        if cds[i][0]<=p<=cds[i][1]: return True
        if p-cds[i][0]>3_000_000: break
        i-=1
    return False
keep=[]; reasons={"maf":0,"lift":0,"cds":0}
for j,k in enumerate(keys):
    if mafT[j]<0.001: reasons["maf"]+=1; continue
    if k not in k2p: reasons["lift"]+=1; continue
    if incds(k2p[k]): reasons["cds"]+=1; continue
    keep.append(j)
n1=len(keep)
if n1>CAP:
    rng=np.random.default_rng(SEED+int(L[1:])); keep=sorted(rng.choice(keep,CAP,replace=False).tolist())
D=D[:,keep].astype(np.float16)
np.savez(f"{W}/private/loci/{L}.npz.tmp.npz", D=D)
os.replace(f"{W}/private/loci/{L}.npz.tmp.npz", f"{W}/private/loci/{L}.npz")
with open(f"{W}/out/loci_meta/{L}.tsv","w") as g:
    g.write("key\tchr\tpos19\tpos38\tmaf_train\tr2\ttyped\trelpos\n")
    for j in keep:
        m=meta[j]; g.write(f"{keys[j]}\t{ch}\t{m[0]}\t{k2p[keys[j]]}\t{mafT[j]:.6g}\t{m[4]}\t{1 if m[5]=='1' else 0}\t{(m[0]-s)/(e-s):.5f}\n")
print(json.dumps({"locus":L,"chr":ch,"span_kb":(e-s)/1e3,"n_raw":n0,"drop":reasons,"n_pass":n1,"n_kept":len(keep)}))
