#!/usr/bin/env python3
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import sys, json, math, time, statistics as st
import numpy as np, scipy.stats as sst
sys.path.insert(0, _config_path("${PROJECT_ROOT}/work/run_l3")); import l3_common as C
C.load_libraries(8)
ROOT=_config_path("${PROJECT_ROOT}/work/ref"); RES=f"{ROOT}/annot/resid"; OUT=_config_path("${PROJECT_ROOT}/work/run_l3b/out")
SEED=20260910; THR=[2.5e-6,1e-5,1e-4,1e-3]; TAUS={"t4":4.0,"t3p5":3.5,"inf":float("inf")}; KS=(2,4); NPERM=int(sys.argv[1]) if len(sys.argv)>1 else 200
t0=time.time()
D=C.Dosage(ROOT,19); samples=D.samples
R=C.Residuals(RES,["tchl"],samples,min_n=1000,shuffle_seed=None)
r_full=R.full[0].astype(np.float64); mask=np.isfinite(r_full)
idx=np.where(mask)[0]; rng=np.random.default_rng(SEED); rng.shuffle(idx)
folds={2:[np.sort(idx[:len(idx)//2]), np.sort(idx[len(idx)//2:])],
       4:[np.sort(f) for f in np.array_split(idx,4)]}
for K in KS:
    fs=folds[K]; u=np.concatenate(fs); assert len(np.unique(u))==len(u)==len(idx), "GATE FAIL folds"
assert np.array_equal(folds[2][0], np.sort(idx[:len(idx)//2])), "GATE FAIL K=2 != ceiling A/B"
comp={K:[np.sort(np.setdiff1d(idx,f)) for f in folds[K]] for K in KS}
for K in KS:
    for f,cmp_ in zip(folds[K],comp[K]): assert len(np.intersect1d(f,cmp_))==0, "GATE FAIL weight/test overlap"
print(f"dosage {D.n_variants} | observed {len(idx)} | K=2 {[len(f) for f in folds[2]]} K=4 {[len(f) for f in folds[4]]} | NPERM {NPERM} | {time.time()-t0:.0f}s",flush=True)
genes={}
for line in open(f"{OUT}/G_phi.txt"):
    t=line.split()
    if len(t)>=3 and t[1]=="var": genes[t[0]]=t[2:]
def center(X): return X-X.mean(axis=1,keepdims=True)
def zfrom(Xc,r):
    rc=r-r.mean(); U=Xc@rc; V=(rc@rc)/len(rc)*np.einsum("ij,ij->i",Xc,Xc)
    return np.where(V>0,U/np.sqrt(np.maximum(V,1e-300)),0.0)
def liu_p(Q,lam):
    lam=lam[lam>1e-10]
    if len(lam)==0: return float("nan")
    c1,c2,c3,c4=[float(np.sum(lam**k)) for k in (1,2,3,4)]; s1=c3/c2**1.5; s2=c4/c2**2
    if s1*s1>s2: a=1/(s1-math.sqrt(s1*s1-s2)); d=s1*a**3-a*a; l=a*a-2*d
    else: a=1/s1; d=0.0; l=1/(s1*s1)
    tstar=(Q-c1)/math.sqrt(2*c2)*math.sqrt(2)*a+l+d
    return float(sst.ncx2.sf(tstar,l,d)) if d>0 else float(sst.chi2.sf(tstar,l))
def eig_w(LD,w): Aw=(w[:,None]*LD)*w[None,:]; return np.linalg.eigvalsh((Aw+Aw.T)/2)
def beta125(maf): return 25.0*(1-maf)**24
def gw(zc,tau,wb): return np.where(np.abs(zc)>=tau,np.abs(zc),1.0)*wb
tracked=["ENSG00000129353.15","ENSG00000213892.12","ENSG00000186567.14","ENSG00000130202.10","ENSG00000130204.13","ENSG00000104856.15","ENSG00000069399.15","ENSG00000079805.19","ENSG00000142453.13","ENSG00000127616.22","ENSG00000129354.12"]
order=[g for g in tracked if g in genes]+[g for g in genes if g not in tracked]
ARMS=["full_flat","full_beta"]+[f"cf{K}_{tn}" for K in KS for tn in TAUS]
PERM_ARMS=["full_beta","cf2_t4","cf4_t4"]
prng=np.random.default_rng(SEED+7); perm_idx=[prng.permutation(len(idx)) for _ in range(NPERM)]
res={}; nullcnt={a:np.zeros((NPERM,len(THR)),dtype=int) for a in PERM_ARMS}
def dump(final):
    cnt={a:[sum(1 for g in res if np.isfinite(res[g][a]) and res[g][a]<t) for t in THR] for a in ARMS}
    def sign(x,y):
        d=[math.log10(res[g][y])-math.log10(res[g][x]) for g in res if all(np.isfinite(res[g][k]) and res[g][k]>0 for k in (x,y))]
        win=sum(1 for v in d if v>1e-12); lose=sum(1 for v in d if v<-1e-12); n=win+lose
        return {"median":st.median(d) if d else None,"win":win,"lose":lose,"tie":len(d)-n,"z":(win-n/2)/math.sqrt(n/4) if n else None}
    pairs=[(a,"full_beta") for a in ARMS if a.startswith("cf")]+[(a,"full_flat") for a in ARMS if a.startswith("cf")]+[("full_beta","full_flat"),("cf2_t4","cf2_inf"),("cf4_t4","cf4_inf")]
    S={f"{x}_vs_{y}":sign(x,y) for x,y in pairs}
    nul={a:{str(t):{"mean":float(nullcnt[a][:,i].mean()),"q95":float(np.quantile(nullcnt[a][:,i],0.95)),"max":int(nullcnt[a][:,i].max()),
             "p_ge_obs":float((np.sum(nullcnt[a][:,i]>=cnt[a][i])+1)/(NPERM+1))} for i,t in enumerate(THR)} for a in PERM_ARMS}
    json.dump({"final":final,"n_genes":len(res),"seed":SEED,"nperm":NPERM,"taus":{k:(v if np.isfinite(v) else "inf") for k,v in TAUS.items()},
               "folds":{str(K):[int(len(f)) for f in folds[K]] for K in KS},"threshold_counts":cnt,"sign_tests":S,"perm_null_counts":nul,"per_gene":res},
              open(f"{OUT}/crossfit{'' if final else '.partial'}.json","w"),indent=1)
    return cnt,S,nul
for gi,g in enumerate(order):
    vloc=D.locate(genes[g]); X=D.M[vloc].toarray().astype(np.float64)
    maf=np.minimum(X.mean(axis=1)/2,1-X.mean(axis=1)/2); wb=beta125(maf); v=len(vloc)
    XcF=center(X[:,idx]); rF=r_full[idx]
    sdF=np.sqrt(np.einsum("ij,ij->i",XcF,XcF)); sdF[sdF==0]=1; LDF=(XcF@XcF.T)/np.outer(sdF,sdF)
    lamF_flat=eig_w(LDF,np.ones(v)); lamF_beta=eig_w(LDF,wb)
    pos={i:j for j,i in enumerate(idx)}
    F={}
    for K in KS:
        for k,(f,cmp_) in enumerate(zip(folds[K],comp[K])):
            cf=np.array([pos[i] for i in f]); cc=np.array([pos[i] for i in cmp_])
            Xk=center(X[:,f]); sd=np.sqrt(np.einsum("ij,ij->i",Xk,Xk)); sd[sd==0]=1
            F[(K,k)]=dict(cf=cf,cc=cc,Xk=Xk,LD=(Xk@Xk.T)/np.outer(sd,sd),Xc=center(X[:,cmp_]))
    def evaluate(r):
        out={}
        zF=zfrom(XcF,r); out["full_flat"]=liu_p(float(np.sum(zF*zF)),lamF_flat); out["full_beta"]=liu_p(float(np.sum(wb*wb*zF*zF)),lamF_beta)
        for K in KS:
            Qs={tn:0.0 for tn in TAUS}; lams={tn:[] for tn in TAUS}
            for k in range(K):
                d=F[(K,k)]; zk=zfrom(d["Xk"],r[d["cf"]]); zc=zfrom(d["Xc"],r[d["cc"]])
                for tn,tau in TAUS.items():
                    w=gw(zc,tau,wb); Qs[tn]+=float(np.sum(w*w*zk*zk)); lams[tn].append(eig_w(d["LD"],w))
            for tn in TAUS: out[f"cf{K}_{tn}"]=liu_p(Qs[tn],np.concatenate(lams[tn]))
        return out
    obs=evaluate(rF); res[g]={"n":int(v),**obs}
    for b in range(NPERM):
        rp=rF[perm_idx[b]]; zF=zfrom(XcF,rp); pb=liu_p(float(np.sum(wb*wb*zF*zF)),lamF_beta)
        pk={}
        for K in KS:
            Q=0.0; lam=[]
            for k in range(K):
                d=F[(K,k)]; zk=zfrom(d["Xk"],rp[d["cf"]]); zc=zfrom(d["Xc"],rp[d["cc"]]); w=gw(zc,4.0,wb)
                Q+=float(np.sum(w*w*zk*zk)); lam.append(eig_w(d["LD"],w))
            pk[K]=liu_p(Q,np.concatenate(lam))
        for i,t in enumerate(THR):
            nullcnt["full_beta"][b,i]+=int(pb<t); nullcnt["cf2_t4"][b,i]+=int(pk[2]<t); nullcnt["cf4_t4"][b,i]+=int(pk[4]<t)
    if g in tracked:
        print(f"  {g:22s} n={v:4d} full_beta={obs['full_beta']:.2e} cf2_t4={obs['cf2_t4']:.2e} cf4_t4={obs['cf4_t4']:.2e} cf2_inf={obs['cf2_inf']:.2e}",flush=True)
    if gi==len(tracked)-1 or (gi%100==0 and gi>0):
        cnt,S,_=dump(False); print(f"[partial {gi+1}/{len(order)} {time.time()-t0:.0f}s] cf2_t4 vs full_beta z={S['cf2_t4_vs_full_beta']['z']}",flush=True)
cnt,S,nul=dump(True)
print("=== 문턱 카운트 (전체 1475 유전자)"); [print(f"  {a:12s}",cnt[a]) for a in ARMS]
print("=== 순열 귀무 카운트 (B=%d)"%NPERM)
for a in PERM_ARMS:
    for t in THR: n_=nul[a][str(t)]; print(f"  {a:10s} <{t:g}: obs {cnt[a][THR.index(t)]} | null mean {n_['mean']:.2f} q95 {n_['q95']:.0f} max {n_['max']} | p_ge {n_['p_ge_obs']:.3f}")
print("=== 부호검정 (동률 제외; 양수 = 앞쪽 우세)")
for k,v_ in S.items(): print(f"  {k:26s} 중위차 {v_['median']:+.4f} | 승 {v_['win']} 패 {v_['lose']} 동률 {v_['tie']} | z {v_['z']:+.2f}")
print(f"CROSSFIT_DONE {time.time()-t0:.0f}s")
