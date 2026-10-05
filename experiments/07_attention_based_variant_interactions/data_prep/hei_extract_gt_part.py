#!/usr/bin/env python
import os as _os
if not _os.environ.get("PROJECT_ROOT", "").strip():
    raise ValueError("PROJECT_ROOT must be set and nonblank")
import sys, os, subprocess, json, numpy as np
N,K=sys.argv[1],int(sys.argv[2]); CH=5000
W=(_os.environ["PROJECT_ROOT"] + '/work/prs'); P=W+"/private/hei"; B=_os.environ.get("BCFTOOLS_BIN", "bcftools")
PP=P+"/gt/parts"; MP=W+"/out/hei/gt/parts"; os.makedirs(PP,exist_ok=True); os.makedirs(MP,exist_ok=True)
tag=f"{N}_{K}"
if os.path.exists(f"{MP}/k_{tag}.done"): sys.exit(0)
keys=[l.strip() for l in open(f"{W}/out/hei/gt/keys_chr{N}.txt") if l.strip()][K*CH:(K+1)*CH]
if not keys: sys.exit(0)
want=set(keys); ids=[l.split("\t")[0] for l in open(P+"/split_hei.tsv")][1:]
sf=f"{PP}/s_{tag}"; open(sf,"w").write("\n".join(ids)+"\n")
rf=f"{PP}/r_{tag}"; open(rf,"w").write("".join(f"{N}\t{k.split(':')[1]}\n" for k in keys))
p=subprocess.Popen([B,"query","-R",rf,"-S",sf,"-f","%POS\t%REF\t%ALT[\t%GT]\n",f"{_os.environ['PROJECT_ROOT']}/work/ref/orig_index/chr{N}.vcf.gz"],stdout=subprocess.PIPE,text=True,bufsize=1<<20)
cols=[]; meta=[]; seen=set()
for line in p.stdout:
    f=line.rstrip("\n").split("\t",3); k=f"{N}:{f[0]}:{f[1]}:{f[2]}"
    if k not in want or k in seen: continue
    seen.add(k); g=f[3].split("\t")
    a=np.array([x[0] for x in g]); b=np.array([x[2] if len(x)>2 else x[0] for x in g])
    unph=sum(1 for x in g if "/" in x)
    h=np.stack([(a=="1"),(b=="1")],axis=1).astype(np.int8); cols.append(h); meta.append((k,unph))
rc=p.wait(); os.remove(sf); os.remove(rf); assert rc==0, rc
H=np.stack(cols,axis=1) if cols else np.zeros((len(ids),0,2),np.int8)
np.save(f"{PP}/h_{tag}.tmp.npy",H); os.replace(f"{PP}/h_{tag}.tmp.npy",f"{PP}/h_{tag}.npy")
with open(f"{MP}/k_{tag}.tsv","w") as o: o.write("key\tn_unphased\n"); [o.write(f"{k}\t{u}\n") for k,u in meta]
json.dump(dict(chr=N,part=K,requested=len(keys),found=len(meta)),open(f"{MP}/k_{tag}.done","w")); print(tag,len(keys),len(meta),flush=True)
