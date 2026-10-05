#!/usr/bin/env python
import os as _os
if not _os.environ.get("PROJECT_ROOT", "").strip():
    raise ValueError("PROJECT_ROOT must be set and nonblank")
import numpy as np, csv, json, collections, os, sys, time
import scipy.linalg as sla
W=(_os.environ["PROJECT_ROOT"] + '/work/prs'); P=W+"/private/hei"; O=W+"/out/hei"
LAMS=np.logspace(-1,9,31); CH=2048
d=np.load(P+"/base_hei.npz",allow_pickle=True); r=d["r"].astype(np.float64); prs=d["prs"].astype(np.float64); sp=d["split"]
tr=sp=="train"; va=sp=="valid"; te=sp=="test"; n=len(r)
A=np.column_stack([np.ones(n),prs]); ba=np.linalg.lstsq(A[tr],r[tr],rcond=None)[0]; base_pred=A@ba; rr=r-base_pred
def r2(y,p): return float(1-((y-p)**2).sum()/((y-y.mean())**2).sum())
loc={}; G={}
for f in sorted(os.listdir(O+"/meta")):
    if not f.startswith("meta_chr"): continue
    N=f[8:-4]; G[N]=np.load(f"{P}/geno_chr{N}.npy",mmap_mode="r")
    for j,row in enumerate(csv.DictReader(open(f"{O}/meta/{f}"),delimiter="\t")): loc[row["key"]]=(N,j)
Ainv=np.linalg.pinv(A[tr])
def std_block(N, cols):
    X=np.asarray(G[N][:,cols],dtype=np.float32)
    mu=X[tr].mean(0); sd=X[tr].std(0); sd[sd<1e-6]=1.0; X=(X-mu)/sd
    B=(Ainv@X[tr].astype(np.float64)).astype(np.float32); X-= (A.astype(np.float32)@B)
    return X
def by_chr(keys):
    m=collections.defaultdict(list)
    for k in keys:
        if k in loc: m[loc[k][0]].append(loc[k][1])
    return {N:sorted(set(v)) for N,v in m.items()}
res={}; pred={"P1":base_pred}
def primal(X,name):
    Xt=X[tr]; Gm=Xt.T@Xt; w,V=np.linalg.eigh(Gm); c=V.T@(Xt.T@rr[tr]); best=None; curve=[]
    for lam in LAMS:
        beta=V@(c/(w+lam)); s=r2(r[va],base_pred[va]+X[va]@beta); curve.append((float(lam),round(s,6)))
        if best is None or s>best[0]: best=(s,lam,beta)
    s,lam,beta=best; res[name]=dict(p=int(X.shape[1]),valid_r2=s,valid_gain=s-r2(r[va],base_pred[va]),**{"lambda":float(lam)},lambda_at_edge=bool(lam in (LAMS[0],LAMS[-1])),curve=curve)
    pred[name]=base_pred+X@beta; print(name,json.dumps({k:v for k,v in res[name].items() if k!="curve"}),flush=True)
def dual(keys,name):
    t0=time.time(); M=by_chr(keys); p=sum(len(v) for v in M.values())
    itr=np.where(tr)[0]; iot=np.where(~tr)[0]
    K=np.zeros((len(itr),len(itr)),np.float64); Ko=np.zeros((len(iot),len(itr)),np.float64)
    for N,cols in M.items():
        for i in range(0,len(cols),CH):
            X=std_block(N,cols[i:i+CH]); Xt=X[itr]
            K+=(Xt@Xt.T); Ko+=(X[iot]@Xt.T); del X,Xt
    print(name,"kernel built p=",p,round(time.time()-t0),"s",flush=True)
    w,U=sla.eigh(K,driver="evr",overwrite_a=True,check_finite=False); del K
    assert w.max()>1.0, f"dual eig degenerate: max={w.max()}"
    c=U.T@rr[itr]; best=None; curve=[]; vmask=va[iot]
    for lam in LAMS:
        a=U@(c/(w+lam)); po=Ko@a; s=r2(r[va],base_pred[va]+po[vmask]); curve.append((float(lam),round(s,6)))
        if best is None or s>best[0]: best=(s,lam,a)
    s,lam,a=best; out=base_pred.copy(); out[iot]+=Ko@a; Kt=None
    res[name]=dict(p=int(p),valid_r2=s,valid_gain=s-r2(r[va],base_pred[va]),**{"lambda":float(lam)},lambda_at_edge=bool(lam in (LAMS[0],LAMS[-1])),eig_min=float(w.min()),eig_max=float(w.max()),curve=curve,sec=round(time.time()-t0))
    pred[name]=out; print(name,json.dumps({k:v for k,v in res[name].items() if k!="curve"}),flush=True)
res["P1"]=dict(valid_r2=r2(r[va],base_pred[va]))
gt=list(csv.DictReader(open(O+"/sel/gene_tokens.tsv"),delimiter="\t")); dt=list(csv.DictReader(open(O+"/sel/dist_tokens.tsv"),delimiter="\t"))
leads=[l.split()[:2] for l in open(O+"/lead_snps.tsv")]
bychr_keys=collections.defaultdict(list)
for k in loc: bychr_keys[k.split(":")[0]].append((int(k.split(":")[1]),k))
for v in bychr_keys.values(): v.sort()
prox=[]
for ch,bp in leads:
    v=bychr_keys.get(ch); 
    if not v: continue
    ps=np.array([x[0] for x in v]); j=int(np.argmin(np.abs(ps-int(bp))))
    if abs(ps[j]-int(bp))<=5000: prox.append(v[j][1])
prox=sorted(set(prox)); M=by_chr(prox)
X2=np.column_stack([std_block(N,cols) for N,cols in M.items()]); res["P2_n_leads"]=len(leads); primal(X2.astype(np.float64),"P2"); del X2
gk=collections.defaultdict(list)
for row in gt: gk[row["gene"]].append(row["key19"])
feats=[]
for g_,ks in gk.items():
    M=by_chr(ks); Xs=[std_block(N,cols) for N,cols in M.items()]
    if not Xs: continue
    X=np.column_stack(Xs); feats.append(X.mean(1))
X4=np.column_stack(feats).astype(np.float64); primal(X4,"P4"); del X4,feats
json.dump(res,open(O+"/train_linear_hei.partial.json","w"))
dual(sorted({row["key19"] for row in gt}),"P3_ANN")
json.dump(res,open(O+"/train_linear_hei.partial.json","w"))
dual(sorted({row["key19"] for row in dt}),"P3_DIST")
np.savez(P+"/pred_linear_hei.npz",**{k:v.astype(np.float32) for k,v in pred.items()})
json.dump(res,open(O+"/train_linear_hei.json.tmp","w"),indent=1); os.replace(O+"/train_linear_hei.json.tmp",O+"/train_linear_hei.json")
print("DONE",flush=True)
