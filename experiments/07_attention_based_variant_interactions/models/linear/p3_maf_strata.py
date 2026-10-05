#!/usr/bin/env python
import sys, json, numpy as np, scipy.linalg as sla
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "data_prep"))
from prs_common import *
TAG=sys.argv[1] if len(sys.argv)>1 else "bbj"
LAMS=np.logspace(-1,9,31)
r,prs,sp=base(TAG); tr=sp=="train"; va=sp=="valid"
A=np.column_stack([np.ones(len(r)),prs]); ba=np.linalg.lstsq(A[tr],r[tr],rcond=None)[0]; bp=A@ba; rr=r-bp
def stratum(f):
    f=min(f,1-f); return "common" if f>=0.05 else "low" if f>=0.01 else "rare"
Ls=loci(); cols={"common":[],"low":[],"rare":[]}
for L in Ls:
    m=meta(L); D=dosage(L,tr).astype(np.float32); s=np.array([stratum(float(x["maf_train"])) for x in m])
    for k in cols: 
        if (s==k).any(): cols[k].append(D[:,s==k])
def fit(X,name):
    X=X.astype(np.float64); B=np.linalg.lstsq(A[tr],X[tr],rcond=None)[0]; X-=A@B
    Xt=X[tr]; p=X.shape[1]; assert p<Xt.shape[0]
    G=Xt.T@Xt; w,V=sla.eigh(G,driver="evr",overwrite_a=True,check_finite=False); del G
    assert w.max()>1.0; c=V.T@(Xt.T@rr[tr]); curve=[]; best=None
    for lam in LAMS:
        beta=V@(c/(w+lam)); s=r2(r[va],bp[va]+X[va]@beta); curve.append((float(lam),round(s,6)))
        if best is None or s>best[0]: best=(s,lam)
    out={"p":int(p),"valid_r2":float(best[0]),"lambda":float(best[1]),"lambda_at_edge":bool(best[1] in (LAMS[0],LAMS[-1])),"eig_max":float(w.max()),"curve":curve}
    print(name, json.dumps({k:v for k,v in out.items() if k!="curve"}), flush=True); return out
res={"tag":TAG,"diagnostic":True,"P1_valid_r2":float(r2(r[va],bp[va]))}
for k in ("rare","low","common"):
    res[k]=fit(np.hstack(cols[k]),k)
res["common_low"]=fit(np.hstack(cols["common"]+cols["low"]),"common_low")
p3=json.load(open(f"{W}/out/train_linear_v3_{TAG}.json"))["P3"]["valid_r2"]; res["P3_all_valid_r2"]=p3
res["rare_increment_over_common_low"]=p3-res["common_low"]["valid_r2"]
tmp=f"{W}/out/p3_maf_strata_{TAG}.json.tmp"; json.dump(res,open(tmp,"w"),indent=1); import os; os.replace(tmp,tmp[:-4])
print("DONE", json.dumps({k:(v["valid_r2"] if isinstance(v,dict) else v) for k,v in res.items()}))
