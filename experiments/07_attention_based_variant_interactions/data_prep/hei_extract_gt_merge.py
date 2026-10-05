#!/usr/bin/env python
import os as _os
if not _os.environ.get("PROJECT_ROOT", "").strip():
    raise ValueError("PROJECT_ROOT must be set and nonblank")
import os, glob, json, numpy as np
W=(_os.environ["PROJECT_ROOT"] + '/work/prs'); P=W+"/private/hei/gt"; PP=P+"/parts"; MP=W+"/out/hei/gt/parts"; tot={}
for kf in sorted(glob.glob(W+"/out/hei/gt/keys_chr*.txt")):
    N=kf.split("keys_chr")[1][:-4]; nk=sum(1 for l in open(kf) if l.strip()); npt=(nk+4999)//5000
    if os.path.exists(f"{W}/out/hei/gt/chr{N}.done"): continue
    assert all(os.path.exists(f"{MP}/k_{N}_{k}.done") for k in range(npt)), f"chr{N} incomplete"
    H=np.concatenate([np.load(f"{PP}/h_{N}_{k}.npy") for k in range(npt)],axis=1)
    np.save(f"{P}/h_chr{N}.tmp.npy",H); os.replace(f"{P}/h_chr{N}.tmp.npy",f"{P}/h_chr{N}.npy")
    L=["key\tn_unphased\n"]
    for k in range(npt): L+=open(f"{MP}/k_{N}_{k}.tsv").readlines()[1:]
    open(f"{W}/out/hei/gt/k_chr{N}.tsv","w").writelines(L)
    unph=sum(int(l.split("\t")[1]) for l in L[1:])
    json.dump(dict(chr=N,requested=nk,found=H.shape[1],unphased_cells=unph),open(f"{W}/out/hei/gt/chr{N}.done","w")); tot[N]=(nk,H.shape[1],unph); del H
    for k in range(npt): os.remove(f"{PP}/h_{N}_{k}.npy")
print(json.dumps(tot))
