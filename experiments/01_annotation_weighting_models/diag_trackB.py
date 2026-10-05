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


import numpy as np, json, math, sys, argparse, time
ap=argparse.ArgumentParser()
ap.add_argument("--scores", default=_config_path("${PROJECT_ROOT}/work/run_trackB/scores_full.tsv"))
ap.add_argument("--variants", default=_config_path("${PROJECT_ROOT}/work/run_trackB/variants_ordered.tsv"))
ap.add_argument("--set", choices=["target","all","ctrl"], default="target")
ap.add_argument("--B", type=int, default=2000)
ap.add_argument("--out", required=True)
ap.add_argument("--min_gene_n", type=int, default=3)
a=ap.parse_args()
t0=time.time()
rng=np.random.default_rng(20260910)
W=_config_path("${PROJECT_ROOT}/work")
S={}; n_bad={"OOB":0,"MISMATCH":0}
with open(a.scores) as fh:
    hdr=fh.readline().rstrip("\n").split("\t")
    for line in fh:
        t=line.rstrip("\n").split("\t")
        if len(t)<len(hdr): continue
        d=dict(zip(hdr,t))
        if d["ref_match"]!="OK": n_bad[d["ref_match"]]=n_bad.get(d["ref_match"],0)+1; continue
        S[d["key"]]=d
pairs=[]
with open(a.variants) as fh:
    hdr=fh.readline().rstrip("\n").split("\t")
    for line in fh:
        d=dict(zip(hdr,line.rstrip("\n").split("\t")))
        if d["key"] not in S: continue
        for g in d["genes"].split(";"):
            pairs.append((d["key"],g,int(d["is_target"]),float(d["z_abs"]),float(d["beta"]),float(d["se"]),int(d["is_indel"])))
import json as _j
vm=_j.load(open(_config_path("${PROJECT_ROOT}/work/run_trackB/variants_manifest.json")))
T=set(vm["target_genes"]); C=set(vm["ctrl_genes"])
sel = T if a.set=="target" else (C if a.set=="ctrl" else T|C)
pairs=[p for p in pairs if p[1] in sel]
per_gene_scored={}
for p in pairs: per_gene_scored[p[1]]=per_gene_scored.get(p[1],0)+1
complete={g:(per_gene_scored.get(g,0), vm["per_gene_n"][g]) for g in sel}
n_complete=sum(1 for g,(s,n) in complete.items() if s>=n*0.98)
key=np.array([p[0] for p in pairs]); gene=np.array([p[1] for p in pairs])
y=np.array([p[3] for p in pairs]); beta=np.array([p[4] for p in pairs]); se=np.array([p[5] for p in pairs]); indel=np.array([p[6] for p in pairs])
zsigned=beta/se
print(json.dumps({"set":a.set,"n_scored_variants":len(S),"n_bad":n_bad,"n_pairs":len(pairs),"n_genes":len(set(gene)),"genes_complete_ge98pct":n_complete,"n_indel_pairs":int(indel.sum())}))
def col(name, f=float): return np.array([f(S[k][name]) for k in key])
M={}
for n in ["S1_all_c3_log","S2_all_full_log","S3_tis_c3_log","S4_tis_full_log","S1r_all_c3_raw","S2r_all_full_raw"]: M[n]=col(n)
M["sign_strength_cage_c3"]=np.abs(col("sign_cage_c3")-0.5)
M["abs_cage_tis_c3_signed"]=np.abs(col("cage_tis_c3_signed"))
M["ref_c3_log"]=col("ref_c3_log")
primary=["S1_all_c3_log","S2_all_full_log","S3_tis_c3_log","S4_tis_full_log","sign_strength_cage_c3"]
z=np.load(f"{W}/ref/annot/cache/fm_all.npz", allow_pickle=True)
m=z["chr"]==19; X=z["X"][m].astype(np.float64); k37=z["key37"][m]; g37=z["gene"][m]; cols=[str(c) for c in z["cols"]]
idx={(str(kk),str(gg)):i for i,(kk,gg) in enumerate(zip(k37,g37))}
rowi=np.array([idx.get((k,g),-1) for k,g in zip(key,gene)])
print(json.dumps({"fm_all_pair_match":int((rowi>=0).sum()),"of":len(rowi)}))
exclude={"maf","r2","avg_cs","is_typed","is_indel","n_genes","gh_n_genes","t1_na"}
cmp_names=[c for c in cols if c not in exclude]
for cname in cmp_names:
    v=np.full(len(key),np.nan); ok=rowi>=0; v[ok]=X[rowi[ok],cols.index(cname)]; M["cmp:"+cname]=v
gidx={}
for i,g in enumerate(gene): gidx.setdefault(g,[]).append(i)
groups=[np.array(v) for v in gidx.values() if len(v)>=a.min_gene_n]
def within_rank(v):
    out=np.full(len(v),np.nan)
    for ix in groups:
        s=v[ix]; ok=~np.isnan(s)
        if ok.sum()>=a.min_gene_n:
            r=np.empty(ok.sum()); r[np.argsort(s[ok],kind="stable")]=np.arange(ok.sum())
            out[ix[ok]]=r-r.mean()
    return out
def corr(x,yy,minn=30):
    ok=~np.isnan(x)&~np.isnan(yy)
    if ok.sum()<minn: return np.nan,int(ok.sum())
    x,yy=x[ok],yy[ok]; x=x-x.mean(); yy=yy-yy.mean()
    den=math.sqrt((x*x).sum()*(yy*yy).sum())
    return (float((x*yy).sum()/den) if den>0 else np.nan), int(ok.sum())
B=a.B
ry=within_rank(y)
ynull=np.empty((B,len(y)))
for b in range(B):
    yy=y.copy()
    for ix in groups: yy[ix]=yy[rng.permutation(ix)]
    ynull[b]=within_rank(yy)
def test(x, ry_, ynull_):
    rx=within_rank(x); rho,nn=corr(rx,ry_)
    if np.isnan(rho): return rho,nn,np.nan,np.nan
    nul=np.array([corr(rx,ynull_[b])[0] for b in range(B)])
    p=(np.sum(np.abs(nul)>=abs(rho))+1)/(B+1)
    return rho,nn,float(np.nanstd(nul)),float(p)
res=[]
for n,v in M.items():
    rho,nn,sd,p=test(v,ry,ynull)
    res.append({"col":n,"rho":None if np.isnan(rho) else rho,"n":nn,"null_sd":None if np.isnan(sd) else sd,"p_perm":None if np.isnan(p) else p,
                "primary":n in primary,"pass_B2":(not np.isnan(rho)) and abs(rho)>=0.05 and p<0.05})
snv=indel==0
groups_all=groups
gidx2={}
for i,g in enumerate(gene):
    if snv[i]: gidx2.setdefault(g,[]).append(i)
groups=[np.array(v) for v in gidx2.values() if len(v)>=a.min_gene_n]
ry_s=within_rank(y); ynull_s=np.empty((B,len(y)))
for b in range(B):
    yy=y.copy()
    for ix in groups: yy[ix]=yy[rng.permutation(ix)]
    ynull_s[b]=within_rank(yy)
supp_snv=[]
for n in primary+["cmp:dist_tss"]:
    rho,nn,sd,p=test(M[n],ry_s,ynull_s); supp_snv.append({"col":n,"rho":None if np.isnan(rho) else rho,"n":nn,"null_sd":None if np.isnan(sd) else sd,"p_perm":None if np.isnan(p) else p})
groups=groups_all
rz=within_rank(zsigned); ysn=np.empty((B,len(y)))
for b in range(B):
    yy=zsigned.copy()
    for ix in groups: yy[ix]=yy[rng.permutation(ix)]
    ysn[b]=within_rank(yy)
sgn_signed=col("cage_tis_c3_signed"); maj=np.sign(col("sign_cage_c3")-0.5)
rho,nn,sd,p=test(sgn_signed,rz,ysn)
supp_sign={"cage_tis_c3_signed_vs_signed_z":{"rho":None if np.isnan(rho) else rho,"n":nn,"null_sd":None if np.isnan(sd) else sd,"p_perm":None if np.isnan(p) else p}}
okm=maj!=0
supp_sign["majority_sign_agree_with_beta_sign"]={"frac":float((maj[okm]==np.sign(beta[okm])).mean()),"n":int(okm.sum())}
okh=okm&(y>=2.0)
supp_sign["majority_sign_agree_with_beta_sign_absz_ge2"]={"frac":float((maj[okh]==np.sign(beta[okh])).mean()) if okh.sum() else None,"n":int(okh.sum())}
def spearman(x,yy):
    ok=~np.isnan(x)&~np.isnan(yy); x,yy=x[ok],yy[ok]
    if len(x)<5: return None
    rx=np.empty(len(x)); rx[np.argsort(x,kind="stable")]=np.arange(len(x)); ryy=np.empty(len(x)); ryy[np.argsort(yy,kind="stable")]=np.arange(len(x))
    r,_=corr(rx,ryy,minn=5); return None if np.isnan(r) else r
per_gene=[]
for g,ix in gidx.items():
    ix=np.array(ix); per_gene.append({"gene":g,"n":len(ix),"is_target":g in T,"rho_S1":spearman(M["S1_all_c3_log"][ix],y[ix]),"rho_S3":spearman(M["S3_tis_c3_log"][ix],y[ix]),"rho_dist_tss":spearman(M["cmp:dist_tss"][ix],y[ix]),"max_absz":float(y[ix].max())})
res.sort(key=lambda t:-abs(t["rho"]) if t["rho"] is not None else 0)
out={"set":a.set,"B":B,"seed":20260910,"rule":"|rho|>=0.05 & p_perm<0.05 (pre-registered B-2)","n_pairs":len(pairs),"n_variants":len(set(key)),"n_genes":len(gidx),"n_groups_ge3":len(groups),
     "genes_complete_ge98pct":n_complete,"per_gene_scored_vs_expected":complete,"n_bad":n_bad,"n_indel_pairs":int(indel.sum()),
     "rows":res,"primary":primary,"passed_primary":[r["col"] for r in res if r["primary"] and r["pass_B2"]],"passed_any":[r["col"] for r in res if r["pass_B2"]],
     "supp_snv_only":supp_snv,"supp_sign":supp_sign,"per_gene":per_gene,"elapsed_s":time.time()-t0}
json.dump(out,open(a.out,"w"),indent=1)
print(f"{'col':28s} {'rho':>8s} {'n':>6s} {'null_sd':>8s} {'p_perm':>7s} prim pass")
for r in res: print(f"{r['col']:28s} {r['rho'] if r['rho'] is not None else float('nan'):+8.4f} {r['n']:6d} {r['null_sd'] if r['null_sd'] is not None else float('nan'):8.4f} {r['p_perm'] if r['p_perm'] is not None else float('nan'):7.4f} {'P' if r['primary'] else ' '}    {'PASS' if r['pass_B2'] else ''}")
print("passed_primary:",out["passed_primary"],"passed_any:",out["passed_any"])
print("supp_sign:",json.dumps(supp_sign)); print("elapsed",round(time.time()-t0,1)); print("DIAG_DONE")
