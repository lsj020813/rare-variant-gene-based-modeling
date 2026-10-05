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


import sys, json, math, csv, statistics as st, time
import numpy as np, scipy.stats as sst
sys.path.insert(0, _config_path("${PROJECT_ROOT}/work/run_l3"))
import l3_common as C
C.load_libraries(8)
ROOT=_config_path("${PROJECT_ROOT}/work/ref"); RES=f"{ROOT}/annot/resid"; OUT=_config_path("${PROJECT_ROOT}/work/run_l3b/out")
SEED=20260910; THR=[2.5e-6,1e-5,1e-4,1e-3]
t0=time.time()
D=C.Dosage(ROOT, 19); samples=D.samples
R=C.Residuals(RES, ["tchl"], samples, min_n=1000, shuffle_seed=None)
r_full=R.full[0].astype(np.float64); mask=np.isfinite(r_full)
idx=np.where(mask)[0]; rng=np.random.default_rng(SEED); rng.shuffle(idx)
A=np.sort(idx[:len(idx)//2]); B=np.sort(idx[len(idx)//2:])
print(f"dosage {D.n_variants} variants, nnz {D.nnz} | tchl observed {len(idx)} | A {len(A)} B {len(B)} | load {time.time()-t0:.0f}s", flush=True)
genes={}
for line in open(f"{OUT}/G_phi.txt"):
    t=line.split()
    if len(t)>=3 and t[1]=="var": genes[t[0]]=t[2:]
def scores(X, r):
    Xc=X-X.mean(axis=1, keepdims=True); rc=r-r.mean()
    U=Xc@rc; V=(rc@rc)/len(rc)*np.einsum("ij,ij->i",Xc,Xc)
    z=np.where(V>0, U/np.sqrt(np.maximum(V,1e-300)), 0.0)
    return z, Xc
def liu_p(Q, lam):
    lam=lam[lam>1e-10]
    if len(lam)==0: return float("nan")
    c1,c2,c3,c4=[float(np.sum(lam**k)) for k in (1,2,3,4)]
    s1=c3/c2**1.5; s2=c4/c2**2
    if s1*s1>s2:
        a=1/(s1-math.sqrt(s1*s1-s2)); d=s1*a**3-a*a; l=a*a-2*d
    else:
        a=1/s1; d=0.0; l=1/(s1*s1)
    muQ=c1; sQ=math.sqrt(2*c2); muX=l+d; sX=math.sqrt(2)*a
    tstar=(Q-muQ)/sQ*sX+muX
    return float(sst.ncx2.sf(tstar, l, d)) if d>0 else float(sst.chi2.sf(tstar, l))
def gene_p(zB, LD, w):
    Q=float(np.sum(w*w*zB*zB)); Aw=(w[:,None]*LD)*w[None,:]
    lam=np.linalg.eigvalsh((Aw+Aw.T)/2)
    return liu_p(Q, lam)
def beta125(maf): return 25.0*(1-maf)**24
tracked=["ENSG00000129353.15","ENSG00000213892.12","ENSG00000186567.14","ENSG00000130202.10","ENSG00000130204.13","ENSG00000104856.15","ENSG00000069399.15","ENSG00000079805.19","ENSG00000142453.13","ENSG00000127616.22","ENSG00000129354.12"]
order=[g for g in tracked if g in genes]+[g for g in genes if g not in tracked]
res={}; arms=["flat","beta","oracle_half","oracle_full","oracle_half_noBeta"]
def dump(final=False):
    cnt={a:[sum(1 for g in res if res[g][a]<t) for t in THR] for a in arms}
    def sign(x,y):
        d=[math.log10(res[g][y])-math.log10(res[g][x]) for g in res if np.isfinite(res[g][x]) and np.isfinite(res[g][y]) and res[g][x]>0 and res[g][y]>0]
        win=sum(1 for v in d if v>1e-12); lose=sum(1 for v in d if v<-1e-12); n=win+lose
        return {"median":st.median(d) if d else None,"win":win,"lose":lose,"tie":len(d)-n,"z":(win-n/2)/math.sqrt(n/4) if n else None}
    S={f"{x}_vs_{y}":sign(x,y) for x,y in [("oracle_half","flat"),("oracle_half","beta"),("oracle_full","flat"),("oracle_full","beta"),("beta","flat"),("oracle_half_noBeta","flat")]}
    json.dump({"final":final,"n_genes":len(res),"nA":int(len(A)),"nB":int(len(B)),"seed":SEED,"threshold_counts":cnt,"sign_tests":S,"per_gene":res},
              open(f"{OUT}/ceiling_halfsplit{'' if final else '.partial'}.json","w"), indent=1)
    return cnt,S
for gi,g in enumerate(order):
    vloc=D.locate(genes[g]); X=D.M[vloc].toarray().astype(np.float64)
    maf=np.minimum(X.mean(axis=1)/2, 1-X.mean(axis=1)/2); wb=beta125(maf)
    zA,_=scores(X[:,A], r_full[A]); zB,XcB=scores(X[:,B], r_full[B]); zF,_=scores(X[:,mask], r_full[mask])
    sd=np.sqrt(np.einsum("ij,ij->i",XcB,XcB)); sd[sd==0]=1
    LD=(XcB@XcB.T)/np.outer(sd,sd)
    res[g]={"n":int(len(vloc)),
            "flat":gene_p(zB,LD,np.ones(len(vloc))),
            "beta":gene_p(zB,LD,wb),
            "oracle_half":gene_p(zB,LD,np.abs(zA)*wb),
            "oracle_half_noBeta":gene_p(zB,LD,np.abs(zA)),
            "oracle_full":gene_p(zB,LD,np.abs(zF)*wb)}
    if g in tracked:
        r_=res[g]; print(f"  {g:22s} n={r_['n']:4d} flat={r_['flat']:.2e} beta={r_['beta']:.2e} oracle_half={r_['oracle_half']:.2e} oracle_full={r_['oracle_full']:.2e} | max|zA| {np.abs(zA).max():.2f} max|zB| {np.abs(zB).max():.2f}", flush=True)
    if gi==len(tracked)-1 or (gi%200==0 and gi>0):
        cnt,S=dump(False); print(f"[partial {gi+1}/{len(order)} {time.time()-t0:.0f}s] oracle_half vs flat z={S['oracle_half_vs_flat']['z']}", flush=True)
cnt,S=dump(True)
print("=== 문턱 카운트 (B 반쪽, n=%d)"%len(B)); [print(f"  {a:20s}", cnt[a]) for a in arms]
print("=== 부호검정 (동률 제외; 양수 = 앞쪽 우세)")
for k,v in S.items(): print(f"  {k:28s} 중위차 {v['median']:+.4f} | 승 {v['win']} 패 {v['lose']} 동률 {v['tie']} | z {v['z']:+.2f}")
print(f"CEILING_DONE {time.time()-t0:.0f}s")
