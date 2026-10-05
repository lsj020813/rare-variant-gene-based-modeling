#!/usr/bin/env python
import sys, json, numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "data_prep"))
from prs_common import *
TAG=sys.argv[1] if len(sys.argv)>1 else "bbj"
LAMS=np.logspace(-2,4,13)
r,prs,sp=base(TAG); tr=sp=="train"; va=sp=="valid"; te=sp=="test"
Ls=loci()
def resid_on(Z, cols):
    A=np.column_stack([np.ones(len(r)),prs]); B=np.linalg.lstsq(A[tr],Z[tr],rcond=None)[0]; return Z-A@B
A=np.column_stack([np.ones(len(r)),prs]); ba=np.linalg.lstsq(A[tr],r[tr],rcond=None)[0]; base_pred=A@ba; rr=r-base_pred
res={"tag":TAG,"n_loci":len(Ls)}
pred={"P0":np.zeros(len(r))+r[tr].mean(),"P1":base_pred}
def ridge_block(X, name):
    X=resid_on(X.astype(np.float64),None)
    Xt=X[tr]; p=X.shape[1]
    best=None
    if p<=Xt.shape[0]:
        G=Xt.T@Xt; w,V=np.linalg.eigh(G); c=V.T@(Xt.T@rr[tr]); del G
        for lam in LAMS*Xt.shape[0]/1e3:
            beta=V@(c/(w+lam)); s=r2(r[va], base_pred[va]+X[va]@beta)
            if best is None or s>best[0]: best=(s,lam,beta)
    else:
        K=Xt@Xt.T; w,U=np.linalg.eigh(K); del K; c=U.T@rr[tr]
        for lam in LAMS*Xt.shape[0]/1e3:
            beta=Xt.T@(U@(c/(w+lam))); s=r2(r[va], base_pred[va]+X[va]@beta)
            if best is None or s>best[0]: best=(s,lam,beta)
    s,lam,beta=best
    res[name]={"p":int(p),"lambda":float(lam),"valid_r2":float(s),"lambda_at_edge":bool(lam in (LAMS[0]*Xt.shape[0]/1e3, LAMS[-1]*Xt.shape[0]/1e3))}
    return base_pred+X@beta
leads=[tuple(l.split()) for l in open(W+"/out/lead_snps.tsv")]
cols2=[]; cols3=[]; cols4=[]
for L in Ls:
    m=meta(L); D=dosage(L,tr); pos=np.array([int(x["pos19"]) for x in m]); ch=m[0]["chr"]
    for c,bp in leads:
        if c==ch and pos.min()-1<=int(bp)<=pos.max()+1:
            cols2.append(D[:,int(np.argmin(np.abs(pos-int(bp))))])
    cols3.append(D.astype(np.float32))
    Dt=D[tr]; U,S,Vt=np.linalg.svd(Dt[:min(20000,len(Dt))],full_matrices=False)
    pcs=D@Vt[:10].T; cols4.append(np.column_stack([D.mean(1),pcs]))
X2=np.column_stack(cols2); pred["P2"]=ridge_block(X2,"P2")
X4=np.column_stack(cols4); pred["P4"]=ridge_block(X4,"P4")
X3=np.column_stack(cols3); del cols3
pred["P3"]=ridge_block(X3,"P3"); del X3
for k in ("P0","P1"): res[k]={"valid_r2":r2(r[va],pred[k][va])}
np.savez(f"{P}/model/pred_linear_{TAG}.npz", **{k:v.astype(np.float32) for k,v in pred.items()})
json.dump(res, open(f"{W}/out/train_linear_{TAG}.json","w"), indent=1); print(json.dumps(res))
