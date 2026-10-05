#!/usr/bin/env python
import os as _os
import math as _number_math
def _required_number(name, cast, positive=False):
    raw = _os.environ.get(name, "")
    if not raw.strip():
        raise ValueError(name + " must be set and nonblank")
    try:
        value = cast(raw)
    except (ValueError, OverflowError):
        raise ValueError(name + " has an invalid numeric value") from None
    if isinstance(value, float) and not _number_math.isfinite(value):
        raise ValueError(name + " must be finite")
    if positive and value <= 0:
        raise ValueError(name + " must be positive")
    return value
if not _os.environ.get("PROJECT_ROOT", "").strip():
    raise ValueError("PROJECT_ROOT must be set and nonblank")
EXPECTED_BASELINE_R2 = _required_number("EXPECTED_BASELINE_R2", float, False)
import numpy as np, csv, json, collections, os, sys, time
import scipy.linalg as sla
W=(_os.environ["PROJECT_ROOT"] + '/work/prs'); P0=W+"/private/hei"; P=W+"/private/hei_seed"; O=W+"/out/hei_seed"; SEL=W+"/out/hei/sel_seed"
LAM=1e6; CH=2048
d=np.load(P0+"/base_hei.npz",allow_pickle=True); r=d["r"].astype(np.float64); prs=d["prs"].astype(np.float64); sp=d["split"]
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

t0=time.time()
gt=[x for x in csv.DictReader(open(SEL+"/gene_tokens.tsv"),delimiter="\t") if int(x["rank"])<128]
M=by_chr(sorted({x["key19"] for x in gt})); p=sum(len(v) for v in M.values()); assert p==123627, p
itr=np.where(tr)[0]; iot=np.where(~tr)[0]
K=np.zeros((len(itr),len(itr)),np.float64); Ko=np.zeros((len(iot),len(itr)),np.float64)
for N,cols in M.items():
    for i in range(0,len(cols),CH):
        X=std_block(N,cols[i:i+CH]); Xt=X[itr]; K+=(Xt@Xt.T); Ko+=(X[iot]@Xt.T); del X,Xt
print("kernel",round(time.time()-t0),"s",flush=True)
y=rr[itr]
Kl=K.copy(); Kl[np.diag_indices_from(Kl)]+=LAM; a=sla.cho_solve(sla.cho_factor(Kl,overwrite_a=True,check_finite=False),y,check_finite=False); del Kl
out=base_pred.copy(); out[iot]+=Ko@a; del Ko
vr=r2(r[va],out[va]); print("full-fit valid_r2",round(vr,6),flush=True)
ref=np.load(P+"/pred_linear_hei_seed.npz")["P3_ANN128"]
assert abs(vr-EXPECTED_BASELINE_R2)<2e-4, vr
maxdiff=float(np.abs(out[~tr]-ref[~tr]).max())
fold=np.random.default_rng(20260928).permutation(len(itr))%5
oofp=np.zeros(len(itr))
for f in range(5):
    fi=np.where(fold!=f)[0]; ho=np.where(fold==f)[0]
    Kf=K[np.ix_(fi,fi)]; Kf[np.diag_indices_from(Kf)]+=LAM
    af=sla.cho_solve(sla.cho_factor(Kf,overwrite_a=True,check_finite=False),y[fi],check_finite=False); del Kf
    oofp[ho]=K[np.ix_(ho,fi)]@af; print("fold",f,round(time.time()-t0),"s",flush=True)
out[itr]=base_pred[itr]+oofp
res=dict(lam=LAM,p=p,valid_r2_fullfit=vr,max_abs_diff_vs_stored=maxdiff,train_oof_r2=r2(r[tr],out[tr]),train_P1_r2=r2(r[tr],base_pred[tr]),
         valid_P1_r2=r2(r[va],base_pred[va]),sec=round(time.time()-t0))
np.save(P+"/p3ann128_oof.tmp.npy",out.astype(np.float64)); os.replace(P+"/p3ann128_oof.tmp.npy",P+"/p3ann128_oof.npy")
json.dump(res,open(O+"/p3ann128_oof.json.tmp","w"),indent=1); os.replace(O+"/p3ann128_oof.json.tmp",O+"/p3ann128_oof.json"); print(json.dumps(res),flush=True)
