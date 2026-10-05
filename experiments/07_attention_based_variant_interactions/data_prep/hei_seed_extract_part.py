#!/usr/bin/env python
import os as _os
if not _os.environ.get("PROJECT_ROOT", "").strip():
    raise ValueError("PROJECT_ROOT must be set and nonblank")
import sys, os, subprocess, json, numpy as np
N,K=sys.argv[1],int(sys.argv[2]); CH=5000
SMOKE=os.environ.get("HEI_SMOKE")=="1"
TAGD="hei_seed_smoke" if SMOKE else "hei_seed"
W=(_os.environ["PROJECT_ROOT"] + '/work/prs'); P0=W+"/private/hei"; P=W+"/private/"+TAGD; B=_os.environ.get("BCFTOOLS_BIN", "bcftools")
PP=P+"/parts"; MP=W+"/out/"+TAGD+"/meta/parts"; os.makedirs(PP,exist_ok=True); os.makedirs(MP,exist_ok=True)
tag=f"{N}_{K}"
if os.path.exists(f"{MP}/m_{tag}.done"): sys.exit(0)
keys=[l.strip() for l in open(f"{W}/out/hei/sel_seed/extract_chr{N}.keys") if l.strip()][K*CH:(K+1)*CH]
if SMOKE: keys=keys[:50]
assert len(keys)==len(set(keys)), "duplicate keys"
if not keys: sys.exit(0)
want=set(keys)
sp=[l.rstrip("\n").split("\t") for l in open(P0+"/split_hei.tsv")][1:]
ids=[a for a,_ in sp]; tr=np.array([b=="train" for _,b in sp])
sf=f"{PP}/s_{tag}"; open(sf,"w").write("\n".join(ids)+"\n")
rf=f"{PP}/r_{tag}"; open(rf,"w").write("".join(f"{N}\t{k.split(':')[1]}\n" for k in keys))
cmd=[B,"query","-R",rf,"-S",sf,"-f","%POS\t%REF\t%ALT\t%INFO/R2\t%INFO/TYPED[\t%DS]\n",f"{_os.environ['PROJECT_ROOT']}/work/ref/orig_index/chr{N}.vcf.gz"]
p=subprocess.Popen(cmd,stdout=subprocess.PIPE,text=True,bufsize=1<<20)
cols=[]; meta=[]; seen=set()
for line in p.stdout:
    f=line.rstrip("\n").split("\t",5); k=f"{N}:{f[0]}:{f[1]}:{f[2]}"
    if k not in want or k in seen: continue
    seen.add(k); d=np.array(f[5].split("\t"),dtype=np.float32); cols.append(d.astype(np.float16))
    af=float(d[tr].mean()/2); meta.append((k,round(min(af,1-af),6),f[3],f[4]))
rc=p.wait(); os.remove(sf); os.remove(rf); assert rc==0, rc
assert len(meta)==len(keys), f"found {len(meta)} of {len(keys)} (missing e.g. {sorted(want-seen)[:3]})"
G=np.stack(cols,axis=1) if cols else np.zeros((len(ids),0),np.float16)
np.save(f"{PP}/g_{tag}.tmp.npy",G); os.replace(f"{PP}/g_{tag}.tmp.npy",f"{PP}/g_{tag}.npy")
with open(f"{MP}/m_{tag}.tsv","w") as g:
    g.write("key\tmaf_train\tr2\ttyped\n"); [g.write("\t".join(map(str,m))+"\n") for m in meta]
json.dump(dict(chr=N,part=K,requested=len(keys),found=len(meta)),open(f"{MP}/m_{tag}.done","w")); print(tag,len(keys),len(meta),flush=True)