#!/usr/bin/env python
import os as _os
if not _os.environ.get("PROJECT_ROOT", "").strip():
    raise ValueError("PROJECT_ROOT must be set and nonblank")
import os, glob, json, numpy as np
SMOKE=os.environ.get("HEI_SMOKE")=="1"
TAGD="hei_seed_smoke" if SMOKE else "hei_seed"
W=(_os.environ["PROJECT_ROOT"] + '/work/prs'); P=W+"/private/"+TAGD+"/gt"; PP=P+"/parts"; MP=W+"/out/"+TAGD+"/gt/parts"; O=W+"/out/"+TAGD+"/gt"; tot={}
for kf in sorted(glob.glob(O+"/keys_chr*.txt")):
    N=kf.split("keys_chr")[1][:-4]; keys=[l.strip() for l in open(kf) if l.strip()]; nk=len(keys); npt=(nk+4999)//5000
    if SMOKE:
        if N!="22": continue
        npt=1; keys=keys[:50]; nk=50
    if os.path.exists(f"{O}/chr{N}.done"): continue
    assert not os.path.exists(f"{P}/h_chr{N}.npy"), f"h_chr{N} exists without done"
    assert all(os.path.exists(f"{MP}/k_{N}_{k}.done") for k in range(npt)), f"chr{N} incomplete"
    H=np.concatenate([np.load(f"{PP}/h_{N}_{k}.npy") for k in range(npt)],axis=1)
    np.save(f"{P}/h_chr{N}.tmp.npy",H); os.replace(f"{P}/h_chr{N}.tmp.npy",f"{P}/h_chr{N}.npy")
    L=["key\tn_unphased\n"]
    for k in range(npt): L+=open(f"{MP}/k_{N}_{k}.tsv").readlines()[1:]
    got=[l.split("\t")[0] for l in L[1:]]
    assert len(got)==H.shape[1] and set(got)==set(keys) and len(set(got))==len(got), f"chr{N} key set mismatch"
    open(f"{O}/k_chr{N}.tsv.tmp","w").writelines(L); os.replace(f"{O}/k_chr{N}.tsv.tmp",f"{O}/k_chr{N}.tsv")
    unph=sum(int(l.split("\t")[1]) for l in L[1:])
    json.dump(dict(chr=N,requested=nk,found=H.shape[1],unphased_cells=unph),open(f"{O}/chr{N}.done","w")); tot[N]=(nk,H.shape[1],unph); del H
    for k in range(npt): os.remove(f"{PP}/h_{N}_{k}.npy")
print(json.dumps(tot))