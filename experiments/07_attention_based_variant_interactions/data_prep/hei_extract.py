#!/usr/bin/env python
import os as _os
if not _os.environ.get("PROJECT_ROOT", "").strip():
    raise ValueError("PROJECT_ROOT must be set and nonblank")
import sys, os, subprocess, json, numpy as np
N=sys.argv[1]; W=(_os.environ["PROJECT_ROOT"] + '/work/prs'); P=W+"/private/hei"; B=_os.environ.get("BCFTOOLS_BIN", "bcftools")
os.makedirs(P,exist_ok=True); os.makedirs(W+"/out/hei/meta",exist_ok=True)
if os.path.exists(f"{W}/out/hei/meta/chr{N}.done"): sys.exit(0)
keys=[l.strip() for l in open(f"{W}/out/hei/sel/extract_chr{N}.keys") if l.strip()]
want=set(keys)
sp=[l.rstrip("\n").split("\t") for l in open(P+"/split_hei.tsv")][1:]
ids=[a for a,_ in sp]; tr=np.array([b=="train" for _,b in sp])
sf=f"{P}/chr{N}.samples"; open(sf,"w").write("\n".join(ids)+"\n")
rf=f"{P}/chr{N}.regions"; open(rf,"w").write("".join(f"{N}\t{k.split(':')[1]}\n" for k in keys))
cmd=[B,"query","-R",rf,"-S",sf,"-f","%POS\t%REF\t%ALT\t%INFO/R2\t%INFO/TYPED[\t%DS]\n",f"{_os.environ['PROJECT_ROOT']}/work/ref/orig_index/chr{N}.vcf.gz"]
p=subprocess.Popen(cmd,stdout=subprocess.PIPE,text=True,bufsize=1<<20)
cols=[]; meta=[]; seen=set()
for line in p.stdout:
    f=line.rstrip("\n").split("\t"); k=f"{N}:{f[0]}:{f[1]}:{f[2]}"
    if k not in want or k in seen: continue
    seen.add(k); d=np.array(f[5:],dtype=np.float32); cols.append(d.astype(np.float16))
    af=float(d[tr].mean()/2); meta.append((k,round(min(af,1-af),6),f[3],f[4]))
rc=p.wait(); os.remove(sf); os.remove(rf)
assert rc==0, f"bcftools rc={rc}"
G=np.stack(cols,axis=1) if cols else np.zeros((len(ids),0),np.float16)
np.save(f"{P}/geno_chr{N}.npy.tmp.npy",G); os.replace(f"{P}/geno_chr{N}.npy.tmp.npy",f"{P}/geno_chr{N}.npy")
with open(f"{W}/out/hei/meta/meta_chr{N}.tsv","w") as g:
    g.write("key\tmaf_train\tr2\ttyped\n"); [g.write("\t".join(map(str,m))+"\n") for m in meta]
json.dump(dict(chr=N,requested=len(keys),found=len(meta),n=len(ids)),open(f"{W}/out/hei/meta/chr{N}.done","w")); print(N,len(keys),len(meta))
