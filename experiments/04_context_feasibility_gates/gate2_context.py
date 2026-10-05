#!/usr/bin/env python
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import gzip, json, os, sys, bisect, collections
import numpy as np
from cyvcf2 import VCF

CHR=sys.argv[1]; GLIST=sys.argv[2]
ROOT=_config_path("${PROJECT_ROOT}/work"); OUT=f"{ROOT}/gate2/out"
ORIG=f"{ROOT}/ref/orig_index/chr{CHR}.vcf.gz"; COMMON=f"{ROOT}/ref/common05/chr{CHR}.maf05.vcf.gz"
CCRE=f"{ROOT}/ref/b6_cards/ccre.s.bed"; PHENO=f"{ROOT}/ref/pheno_v3/tchl_v3.tsv"
G1=f"{ROOT}/gate1/out"
SEED=20260921; NFOLD=5; NPERM=200; MAXPAIR=50; LAMS=[1e-2,1e-1,1.0,10.0,100.0,1000.0]
MIN_COCARRIER=50
rng=np.random.default_rng(SEED)
os.makedirs(OUT, exist_ok=True)

vcf_ids=VCF(ORIG).samples
idx_of={s:i for i,s in enumerate(vcf_ids)}
ph_rows=[]
with open(PHENO) as f:
    hdr=f.readline().rstrip("\n").split("\t")
    ci={k:i for i,k in enumerate(hdr)}
    for line in f:
        p=line.rstrip("\n").split("\t")
        s=p[ci["sample_id"]]
        if s in idx_of: ph_rows.append((idx_of[s], p))
ph_rows.sort()
keep=np.array([r[0] for r in ph_rows])
COVS=["age","sex_male","CT","NC","PC1","PC2","PC3","PC4","PC5"]
y=np.array([float(r[1][ci["y"]]) for r in ph_rows], dtype=np.float64)
Cv=np.column_stack([np.ones(len(ph_rows))]+[np.array([float(r[1][ci[c]]) for r in ph_rows]) for c in COVS])
N=len(y)
fold=rng.permutation(N)%NFOLD

yres=np.zeros(N)
for f_ in range(NFOLD):
    tr=fold!=f_; te=~tr
    beta,*_=np.linalg.lstsq(Cv[tr], y[tr], rcond=None)
    yres[te]=y[te]-Cv[te]@beta
yres-= yres.mean(); SST=float((yres**2).sum())

with gzip.open(f"{G1}/windows_chr{CHR}.json.gz","rt") as f: WIN=json.load(f)
VP={}
with gzip.open(f"{G1}/vpos_chr{CHR}.tsv.gz","rt") as f:
    f.readline()
    for line in f:
        a,b,c_=line.rstrip("\n").split("\t"); VP[a]=(int(b),int(c_))
cc=[]
with open(CCRE) as f:
    for line in f:
        p=line.rstrip("\n").split("\t")
        if p[0] not in (f"chr{CHR}",CHR): continue
        cc.append((int(p[1]),int(p[2])))
cc.sort(); ccs=[x[0] for x in cc]
def ccre_of(p):
    i=bisect.bisect_right(ccs,p)-1
    return i if (i>=0 and cc[i][1]>=p) else -1

def read_dosage(ids):
    recs=sorted((int(v.split(":")[1]), v) for v in ids)
    want=collections.defaultdict(list)
    for pos,v in recs:
        c_,p_,ref,alt=v.split(":"); want[int(p_)].append((ref,alt,v))
    regions=[]
    for pos in sorted(want):
        if regions and pos-regions[-1][1]<100000: regions[-1][1]=pos
        else: regions.append([pos,pos])
    out={}
    for src in (ORIG,COMMON):
        if not os.path.exists(src): continue
        v=VCF(src)
        for s,e in regions:
            try: it=v(f"{CHR}:{max(1,s)}-{e}")
            except Exception: continue
            for rec in it:
                cand=want.get(rec.POS)
                if not cand: continue
                for ref,alt,vid in cand:
                    if vid in out or rec.REF!=ref or alt not in (rec.ALT or []): continue
                    ds=rec.format('DS')
                    if ds is None: continue
                    out[vid]=np.asarray(ds,dtype=np.float32).ravel()[keep]
        v.close()
    kept=[v for _,v in recs if v in out]
    if not kept: return np.zeros((0,N),dtype=np.float32),[]
    X=np.vstack([out[v] for v in kept]); np.nan_to_num(X,copy=False,nan=0.0)
    return X, kept

def ld_prune(X, thr=0.8):
    if X.shape[0]<2: return list(range(X.shape[0]))
    sd=X.std(axis=1); ok=np.where(sd>0)[0]
    if ok.size<2: return list(ok)
    Z=(X[ok]-X[ok].mean(axis=1,keepdims=True))/sd[ok][:,None]
    R2=((Z@Z.T)/X.shape[1])**2
    sel=[]
    for i in range(ok.size):
        if all(R2[i,j]<thr for j in sel): sel.append(i)
    return [int(ok[i]) for i in sel]

def fit_eval(F, lam_grid=LAMS):
    if F.shape[1]==0:
        return SST, [], [(np.zeros((0,0)), np.arange(0))]
    pred=np.zeros(N); ops=[]; lams=[]
    for f_ in range(NFOLD):
        tr=np.where(fold!=f_)[0]; te=np.where(fold==f_)[0]
        Ft=F[tr]; yt=yres[tr]
        G=Ft.T@Ft; b=Ft.T@yt
        inner=rng.permutation(tr.size)%NFOLD
        best=(np.inf,lam_grid[0])
        for lam in lam_grid:
            err=0.0
            for k in range(NFOLD):
                itr=tr[inner!=k]; ite=tr[inner==k]
                Gi=F[itr].T@F[itr]+lam*np.eye(F.shape[1]); bi=F[itr].T@yres[itr]
                try: w=np.linalg.solve(Gi,bi)
                except np.linalg.LinAlgError: w=np.linalg.lstsq(Gi,bi,rcond=None)[0]
                err+=float(((yres[ite]-F[ite]@w)**2).sum())
            if err<best[0]: best=(err,lam)
        lam=best[1]; lams.append(lam)
        A=np.linalg.solve(G+lam*np.eye(F.shape[1]), F[tr].T)
        ops.append((F[te]@A, tr, te))
        pred[te]=F[te]@(A@yt)
    return float(((yres-pred)**2).sum()), lams, ops

def sse_perm(ops, yp):
    pred=np.zeros(N)
    for P,tr,te in ops: pred[te]=P@yp[tr]
    return float(((yp-pred)**2).sum())

genes=[l.strip() for l in open(GLIST) if l.strip()]
res=[]
for gi,g in enumerate(genes):
    ent=WIN.get(g)
    if ent is None: continue
    low=[r[0] for r in ent["bins"].get("B2",[]) if VP.get(r[0],(0,1))[1]==0]
    com=[r[0] for r in ent["bins"].get("B4",[]) if VP.get(r[0],(0,1))[1]==0]
    if len(low)<2 or len(com)<1:
        res.append({"gene":g,"status":"insufficient_variants","n_low":len(low),"n_common":len(com)}); continue
    X,kept=read_dosage(low+com)
    if X.shape[0]<3:
        res.append({"gene":g,"status":"dosage_read_failed"}); continue
    ki={v:i for i,v in enumerate(kept)}
    li=[ki[v] for v in low if v in ki]; cio=[ki[v] for v in com if v in ki]
    if len(li)<2 or len(cio)<1:
        res.append({"gene":g,"status":"insufficient_after_read","n_low":len(li),"n_common":len(cio)}); continue
    Xl=X[li]; Xc=X[cio]
    cpr=ld_prune(Xc, 0.8); Xcp=Xc[cpr]
    flat=Xl.sum(axis=0)[:,None]
    F0=np.column_stack([Xcp.T.astype(np.float64), flat.astype(np.float64)])
    F0=(F0-F0.mean(axis=0))/(F0.std(axis=0)+1e-9)
    sse0,lam0,ops0=fit_eval(F0)
    units=collections.defaultdict(lambda: {"low":[], "com":[]})
    for j,v in enumerate(kept):
        u=ccre_of(VP.get(v,(0,0))[0])
        if u<0: continue
        units[u]["low" if j in set(li) else "com"].append(j)
    ul=[u for u,d in units.items() if d["low"]]
    rec={"gene":g,"chr":CHR,"status":"ok","n_low":len(li),"n_common":len(cio),"n_common_pruned":len(cpr),
         "n_units_with_low":len(ul),"n_units_total":len(units),"SST":SST,"N":N,
         "sse0":sse0,"lambda0":lam0}
    def delta(F, tag):
        Fz=np.column_stack([F0, (F-F.mean(axis=0))/(F.std(axis=0)+1e-9)])
        sse,lam,ops=fit_eval(Fz)
        d=(sse0-sse)/SST
        null=np.empty(NPERM)
        for b in range(NPERM):
            yp=yres[rng.permutation(N)]
            null[b]=((sse_perm(ops0,yp)-sse_perm(ops,yp))/SST)
        rec[f"dR2_{tag}"]=float(d); rec[f"p_perm_{tag}"]=float((null>=d).mean())
        rec[f"null_mean_{tag}"]=float(null.mean()); rec[f"null_sd_{tag}"]=float(null.std())
        rec[f"null_p95_{tag}"]=float(np.quantile(null,0.95)); rec[f"lambda_{tag}"]=lam
        return d
    if ul:
        E=np.vstack([Xl[[li.index(j) for j in units[u]["low"]]].sum(axis=0) for u in ul]).T.astype(np.float64)
        delta(E,"M1")
        sizes=[len(units[u]["low"]) for u in ul]; perm=rng.permutation(len(li)); pos=0; Esh=[]
        for s_ in sizes:
            Esh.append(Xl[perm[pos:pos+s_]].sum(axis=0)); pos+=s_
        Esh=np.vstack(Esh).T.astype(np.float64)
        Fz=np.column_stack([F0,(Esh-Esh.mean(axis=0))/(Esh.std(axis=0)+1e-9)])
        s_sh,_,_=fit_eval(Fz); rec["dR2_NC3_map_shuffle"]=float((sse0-s_sh)/SST)
        Esc=E[rng.permutation(N)]
        Fz=np.column_stack([F0,(Esc-Esc.mean(axis=0))/(Esc.std(axis=0)+1e-9)])
        s_sc,_,_=fit_eval(Fz); rec["dR2_NC2_ctx_shuffle"]=float((sse0-s_sc)/SST)
    if ul:
        inter=[]
        for u in ul:
            if not units[u]["com"]: continue
            lo=Xl[[li.index(j) for j in units[u]["low"]]].sum(axis=0)
            co_=X[units[u]["com"]].sum(axis=0)
            inter.append(lo*(co_-co_.mean()))
        if inter:
            I=np.vstack(inter).T.astype(np.float64); rec["n_inter_units"]=I.shape[1]
            delta(I,"M3")
    Bl=(Xl>=0.5).astype(np.float32); Bc=(Xc>=0.5).astype(np.float32)
    COv=Bl@Bc.T
    cand=np.argwhere(COv>=MIN_COCARRIER)
    if cand.size:
        order=np.argsort(-COv[cand[:,0],cand[:,1]])[:MAXPAIR]
        sel=cand[order]
        P=np.vstack([Xl[a]*Xc[b] for a,b in sel]).T.astype(np.float64)
        rec["n_pairs_selected"]=int(sel.shape[0]); rec["n_pairs_ge50"]=int(cand.shape[0])
        delta(P,"M2")
        pos_l=np.array([int(kept[li[a]].split(":")[1]) for a,_ in sel])
        pos_c=np.array([int(kept[cio[b]].split(":")[1]) for _,b in sel])
        dists=np.abs(pos_l-pos_c)
        allpos_l=np.array([int(kept[j].split(":")[1]) for j in li])
        allpos_c=np.array([int(kept[j].split(":")[1]) for j in cio])
        Dmat=np.abs(allpos_l[:,None]-allpos_c[None,:])
        pick=[]
        for d_ in dists:
            j=np.argmin(np.abs(Dmat-d_)+rng.normal(0,1,Dmat.shape))
            pick.append(np.unravel_index(j,Dmat.shape))
        P4=np.vstack([Xl[a]*Xc[b] for a,b in pick]).T.astype(np.float64)
        Fz=np.column_stack([F0,(P4-P4.mean(axis=0))/(P4.std(axis=0)+1e-9)])
        s4,_,_=fit_eval(Fz); rec["dR2_NC4_dist_matched"]=float((sse0-s4)/SST)
        sdl=Xl.std(axis=1); sdc=Xc.std(axis=1)
        Zl=(Xl-Xl.mean(axis=1,keepdims=True))/(sdl[:,None]+1e-9); Zc=(Xc-Xc.mean(axis=1,keepdims=True))/(sdc[:,None]+1e-9)
        R2lc=((Zl@Zc.T)/Xl.shape[1])**2
        tgt=R2lc[sel[:,0],sel[:,1]]
        pick5=[]
        for t_ in tgt:
            j=np.argmin(np.abs(R2lc-t_)+rng.normal(0,1e-3,R2lc.shape))
            pick5.append(np.unravel_index(j,R2lc.shape))
        P5=np.vstack([Xl[a]*Xc[b] for a,b in pick5]).T.astype(np.float64)
        Fz=np.column_stack([F0,(P5-P5.mean(axis=0))/(P5.std(axis=0)+1e-9)])
        s5,_,_=fit_eval(Fz); rec["dR2_NC5_ld_matched"]=float((sse0-s5)/SST)
        rec["median_cocarrier_selected"]=float(np.median(COv[sel[:,0],sel[:,1]]))
        rec["median_r2_selected"]=float(np.median(tgt))
    res.append(rec)
    if (gi+1)%5==0: print(f"progress chr{CHR} {gi+1}/{len(genes)}", flush=True)

tmp=f"{OUT}/g2_chr{CHR}.json.tmp"
json.dump({"chr":CHR,"seed":SEED,"nfold":NFOLD,"nperm":NPERM,"maxpair":MAXPAIR,
           "min_cocarrier":MIN_COCARRIER,"lambda_grid":LAMS,"N":int(N),"SST":SST,
           "primary_metric":"pooled out-of-fold Delta R^2 vs M0","genes":res}, open(tmp,"w"))
os.replace(tmp,f"{OUT}/g2_chr{CHR}.json")
print(json.dumps({"chr":CHR,"ok":sum(1 for r in res if r.get("status")=="ok"),
                  "failed":sum(1 for r in res if r.get("status")!="ok")}))
