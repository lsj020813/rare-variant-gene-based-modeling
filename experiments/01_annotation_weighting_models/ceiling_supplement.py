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


import sys, json, math, statistics as st, time
sys.path.insert(0, _config_path("${PROJECT_ROOT}/work/run_l3"))
import l3_common as C
C.load_libraries(8)
import numpy as np, scipy.stats as sst
ROOT=_config_path("${PROJECT_ROOT}/work/ref"); RES=f"{ROOT}/annot/resid"; OUT=_config_path("${PROJECT_ROOT}/work/run_l3b/out")
SEED=20260910; THR=[2.5e-6,1e-5,1e-4,1e-3]; NPERM=20
t0=time.time()
D=C.Dosage(ROOT,19); samples=D.samples
R=C.Residuals(RES,["tchl"],samples,min_n=1000,shuffle_seed=None)
r_full=R.full[0].astype(np.float64); mask=np.isfinite(r_full)
idx=np.where(mask)[0]; rng=np.random.default_rng(SEED); rng.shuffle(idx)
A=np.sort(idx[:len(idx)//2]); B=np.sort(idx[len(idx)//2:])
genes={}
for line in open(f"{OUT}/G_phi.txt"):
    t=line.split()
    if len(t)>=3 and t[1]=="var": genes[t[0]]=t[2:]
def scores(X,r):
    Xc=X-X.mean(axis=1,keepdims=True); rc=r-r.mean(); U=Xc@rc; V=(rc@rc)/len(rc)*np.einsum("ij,ij->i",Xc,Xc)
    return np.where(V>0,U/np.sqrt(np.maximum(V,1e-300)),0.0), Xc
def liu_p(Q,lam):
    lam=lam[lam>1e-10]
    if len(lam)==0: return float("nan")
    c1,c2,c3,c4=[float(np.sum(lam**k)) for k in (1,2,3,4)]; s1=c3/c2**1.5; s2=c4/c2**2
    if s1*s1>s2: a=1/(s1-math.sqrt(s1*s1-s2)); d=s1*a**3-a*a; l=a*a-2*d
    else: a=1/s1; d=0.0; l=1/(s1*s1)
    tstar=(Q-c1)/math.sqrt(2*c2)*math.sqrt(2)*a+l+d
    return float(sst.ncx2.sf(tstar,l,d)) if d>0 else float(sst.chi2.sf(tstar,l))
def gene_p(zB,LD,w):
    Q=float(np.sum(w*w*zB*zB)); Aw=(w[:,None]*LD)*w[None,:]
    return liu_p(Q,np.linalg.eigvalsh((Aw+Aw.T)/2))
def beta125(m): return 25.0*(1-m)**24
prng=np.random.default_rng(SEED+1)
res={}
for gi,(g,vk) in enumerate(genes.items()):
    vloc=D.locate(vk); X=D.M[vloc].toarray().astype(np.float64)
    maf=np.minimum(X.mean(axis=1)/2,1-X.mean(axis=1)/2); wb=beta125(maf)
    zA,_=scores(X[:,A],r_full[A]); zB,XcB=scores(X[:,B],r_full[B])
    sd=np.sqrt(np.einsum("ij,ij->i",XcB,XcB)); sd[sd==0]=1; LD=(XcB@XcB.T)/np.outer(sd,sd)
    flat=gene_p(zB,LD,np.ones(len(vloc))); orc=gene_p(zB,LD,np.abs(zA)*wb)
    perms=[]
    for _ in range(NPERM):
        perms.append(gene_p(zB,LD,prng.permutation(np.abs(zA))*wb))
    res[g]={"n":int(len(vloc)),"maxzA":float(np.abs(zA).max()),"flat":flat,"oracle_half":orc,"perm":perms}
    if gi%300==0: print(f"[{gi}/{len(genes)} {time.time()-t0:.0f}s]",flush=True)
def sign(gs,x,y):
    d=[math.log10(res[g][y])-math.log10(res[g][x]) for g in gs if res[g][x]>0 and res[g][y]>0]
    win=sum(1 for v in d if v>1e-12); lose=sum(1 for v in d if v<-1e-12); n=win+lose
    return {"n_genes":len(gs),"median":st.median(d) if d else None,"win":win,"lose":lose,"tie":len(d)-n,"z":(win-n/2)/math.sqrt(n/4) if n>0 else None}
out={"seed":SEED,"nperm":NPERM}
print("=== (1) A 선택 유전자 (B 미열람) — oracle_half vs flat")
for cut in (3.5,4.0,4.5,5.0):
    gs=[g for g in res if res[g]["maxzA"]>=cut]; s=sign(gs,"oracle_half","flat"); out[f"sel_maxzA_ge_{cut}"]=s
    print(f"  max|zA|>={cut}: 유전자 {len(gs):4d} | 중위차 {s['median']:+.3f} | 승 {s['win']} 패 {s['lose']} | z {s['z']:+.2f}" if s['z'] is not None else f"  >={cut}: {len(gs)}")
gs=[g for g in res if res[g]["maxzA"]<3.5]; s=sign(gs,"oracle_half","flat"); out["sel_maxzA_lt_3.5"]=s
print(f"  max|zA|<3.5 (잡음군): 유전자 {len(gs)} | 중위차 {s['median']:+.3f} | 승 {s['win']} 패 {s['lose']} | z {s['z']:+.2f}")
print("=== (2) 문턱 카운트 순열 귀무 (|zA| 유전자 내 셔플 x%d)"%NPERM)
cnt={}
for t in THR:
    obs=sum(1 for g in res if res[g]["oracle_half"]<t); fl=sum(1 for g in res if res[g]["flat"]<t)
    nul=[sum(1 for g in res if res[g]["perm"][k]<t) for k in range(NPERM)]
    cnt[str(t)]={"flat":fl,"oracle_half":obs,"perm_mean":st.mean(nul),"perm_max":max(nul),"perm_min":min(nul),"p_ge":(sum(1 for x in nul if x>=obs)+1)/(NPERM+1)}
    print(f"  <{t:.0e}: flat {fl} | oracle_half {obs} | 순열 평균 {st.mean(nul):.1f} 범위 [{min(nul)},{max(nul)}] | P(순열>=obs) {(sum(1 for x in nul if x>=obs)+1)/(NPERM+1):.3f}")
out["counts"]=cnt
json.dump(out,open(f"{OUT}/ceiling_supplement.json","w"),indent=1)
print(f"CEILSUP_DONE {time.time()-t0:.0f}s")
