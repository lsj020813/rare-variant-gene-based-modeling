#!/usr/bin/env python
import os as _os
if not _os.environ.get("PROJECT_ROOT", "").strip():
    raise ValueError("PROJECT_ROOT must be set and nonblank")
import os, glob, json, numpy as np
W=(_os.environ["PROJECT_ROOT"] + '/work/prs'); P=W+"/private/hei"; PP=P+"/parts"; MP=W+"/out/hei/meta/parts"; tot={}
for kf in sorted(glob.glob(W+"/out/hei/sel/extract_chr*.keys")):
    N=kf.split("extract_chr")[1][:-5]; nk=sum(1 for l in open(kf) if l.strip()); np_=(nk+4999)//5000
    if os.path.exists(f"{W}/out/hei/meta/chr{N}.done"): continue
    assert all(os.path.exists(f"{MP}/m_{N}_{k}.done") for k in range(np_)), f"chr{N} incomplete"
    G=np.concatenate([np.load(f"{PP}/g_{N}_{k}.npy") for k in range(np_)],axis=1)
    np.save(f"{P}/geno_chr{N}.tmp.npy",G); os.replace(f"{P}/geno_chr{N}.tmp.npy",f"{P}/geno_chr{N}.npy")
    lines=["key\tmaf_train\tr2\ttyped\n"]
    for k in range(np_): lines+=open(f"{MP}/m_{N}_{k}.tsv").readlines()[1:]
    open(f"{W}/out/hei/meta/meta_chr{N}.tsv","w").writelines(lines)
    json.dump(dict(chr=N,requested=nk,found=G.shape[1]),open(f"{W}/out/hei/meta/chr{N}.done","w")); tot[N]=(nk,G.shape[1]); del G
    for k in range(np_): os.remove(f"{PP}/g_{N}_{k}.npy")
print(json.dumps(tot))
