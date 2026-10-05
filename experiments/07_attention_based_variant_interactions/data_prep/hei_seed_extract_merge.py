#!/usr/bin/env python
import os as _os
if not _os.environ.get("PROJECT_ROOT", "").strip():
    raise ValueError("PROJECT_ROOT must be set and nonblank")
import os, glob, json, numpy as np
SMOKE=os.environ.get("HEI_SMOKE")=="1"
TAGD="hei_seed_smoke" if SMOKE else "hei_seed"
W=(_os.environ["PROJECT_ROOT"] + '/work/prs'); P=W+"/private/"+TAGD; PP=P+"/parts"; MP=W+"/out/"+TAGD+"/meta/parts"; O=W+"/out/"+TAGD; tot={}
for kf in sorted(glob.glob(W+"/out/hei/sel_seed/extract_chr*.keys")):
    N=kf.split("extract_chr")[1][:-5]; keys=[l.strip() for l in open(kf) if l.strip()]; nk=len(keys); np_=(nk+4999)//5000
    if SMOKE:
        if N!="22": continue
        np_=1; keys=keys[:50]; nk=50
    if os.path.exists(f"{O}/meta/chr{N}.done"): continue
    assert not os.path.exists(f"{P}/geno_chr{N}.npy"), f"geno_chr{N} exists without done; manual review"
    assert all(os.path.exists(f"{MP}/m_{N}_{k}.done") for k in range(np_)), f"chr{N} incomplete"
    G=np.concatenate([np.load(f"{PP}/g_{N}_{k}.npy") for k in range(np_)],axis=1)
    np.save(f"{P}/geno_chr{N}.tmp.npy",G); os.replace(f"{P}/geno_chr{N}.tmp.npy",f"{P}/geno_chr{N}.npy")
    lines=["key\tmaf_train\tr2\ttyped\n"]
    for k in range(np_): lines+=open(f"{MP}/m_{N}_{k}.tsv").readlines()[1:]
    got=[l.split("\t")[0] for l in lines[1:]]
    assert len(got)==G.shape[1] and len(set(got))==len(got) and set(got)==set(keys), f"chr{N} key set mismatch"
    open(f"{O}/meta/meta_chr{N}.tsv.tmp","w").writelines(lines); os.replace(f"{O}/meta/meta_chr{N}.tsv.tmp",f"{O}/meta/meta_chr{N}.tsv")
    json.dump(dict(chr=N,requested=nk,found=G.shape[1]),open(f"{O}/meta/chr{N}.done","w")); tot[N]=(nk,G.shape[1]); del G
    for k in range(np_): os.remove(f"{PP}/g_{N}_{k}.npy")
print(json.dumps(tot))