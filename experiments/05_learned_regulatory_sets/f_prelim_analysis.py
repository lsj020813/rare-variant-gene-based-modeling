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


import gzip, glob, json, os, sys
import numpy as np, pandas as pd
from scipy import stats
from sklearn.cluster import AgglomerativeClustering
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score
from sklearn.neighbors import NearestNeighbors

ROOT=_config_path("${PROJECT_ROOT}/work"); OUT=f"{ROOT}/fset/prelim"; FO=f"{ROOT}/fset/out"
os.makedirs(OUT, exist_ok=True)
SEED=20260922; rng=np.random.default_rng(SEED)

fs=sorted(glob.glob(f"{FO}/feat_chr*.tsv.gz"))
df=pd.concat([pd.read_csv(f, sep="\t", low_memory=False) for f in fs], ignore_index=True)
smp=pd.read_csv(f"{OUT}/prelim_sample.tsv", sep="\t")
df=df.merge(smp[["key","maf","r2","typed","domain_n_total"]], on="key", how="left")

vcols=[c for c in df.columns if c.startswith("re2g")]
cat=["ccre_class","rep_class"]
num=["ccre_dist","rep_dist","cpg","cpg_dist","map_k36","tf_n"]+[c for c in df.columns if c.startswith("tf_") and c!="tf_n"]
X=pd.get_dummies(df[cat].fillna("none"), columns=cat, dtype=float)
for c in num:
    v=pd.to_numeric(df[c], errors="coerce")
    X[f"miss_{c}"]=v.isna().astype(float)
    if c.endswith("_dist"):
        v=np.log10(v.clip(lower=0)+1)
    X[c]=v.fillna(v.median())
keep=[c for c in X.columns if X[c].std()>0]
X=X[keep].astype(np.float32)
qc=pd.DataFrame({"feature":X.columns,
                 "sd":X.std().values,
                 "n_unique":[X[c].nunique() for c in X.columns],
                 "frac_zero":(X==0).mean().values})
qc.to_csv(f"{OUT}/prelim_feature_qc.csv", index=False)
Z=(X-X.mean())/X.std()
Z=Z.values.astype(np.float32)

def hopkins(A, m=500, draws=20, seed=SEED):
    out=[]
    r=np.random.default_rng(seed)
    lo, hi = A.min(axis=0), A.max(axis=0)
    nn=NearestNeighbors(n_neighbors=2).fit(A)
    for _ in range(draws):
        idx=r.choice(A.shape[0], m, replace=False)
        u=r.uniform(lo, hi, size=(m, A.shape[1])).astype(np.float32)
        du=nn.kneighbors(u, n_neighbors=1, return_distance=True)[0][:,0]
        dw=nn.kneighbors(A[idx], n_neighbors=2, return_distance=True)[0][:,1]
        out.append(du.sum()/(du.sum()+dw.sum()))
    return np.array(out)

def eff_rank(A):
    p=PCA(n_components=min(30, A.shape[1])).fit(A)
    ev=p.explained_variance_ratio_
    ent=-(ev[ev>0]*np.log(ev[ev>0])).sum()
    return float(np.exp(ent)), ev

def nn1(A, m=5000, seed=SEED):
    r=np.random.default_rng(seed)
    idx=r.choice(A.shape[0], min(m, A.shape[0]), replace=False)
    nn=NearestNeighbors(n_neighbors=2).fit(A)
    d=nn.kneighbors(A[idx], n_neighbors=2, return_distance=True)[0][:,1]
    return d

H_real=hopkins(Z); er_real, ev_real=eff_rank(Z); d_real=nn1(Z)
H_null=[]; er_null=[]; d_null=[]
for k in range(20):
    S=Z.copy()
    for j in range(S.shape[1]):
        S[:,j]=S[rng.permutation(S.shape[0]),j]
    H_null.append(hopkins(S, draws=3, seed=SEED+k).mean())
    if k<5:
        er_null.append(eff_rank(S)[0]); d_null.append(np.median(nn1(S, seed=SEED+k)))
H_null=np.array(H_null)

def ci(a, q=(2.5,97.5)): return [float(np.percentile(a,q[0])), float(np.percentile(a,q[1]))]
f1={"n_variants":int(Z.shape[0]), "n_features":int(Z.shape[1]),
    "hopkins_real_median":float(np.median(H_real)), "hopkins_real_ci":ci(H_real),
    "hopkins_shuffled_median":float(np.median(H_null)), "hopkins_shuffled_ci":ci(H_null),
    "eff_rank_real":er_real, "eff_rank_shuffled_mean":float(np.mean(er_null)),
    "pc1_var":float(ev_real[0]), "pc_var_top5":[float(x) for x in ev_real[:5]],
    "nn1_median_real":float(np.median(d_real)), "nn1_median_shuffled":float(np.mean(d_null))}

rows=[]
doms=df.domain.value_counts()
for g in doms[doms>=40].index:
    m=df.domain==g
    A=Z[m.values]; pos=df.loc[m,"pos38"].values.astype(float); mafv=df.loc[m,"maf"].values
    if A.shape[0]<40: continue
    n=A.shape[0]
    ii=rng.choice(n, min(n,150), replace=False)
    Asub=A[ii]; psub=pos[ii]
    fd=np.sqrt(((Asub[:,None,:]-Asub[None,:,:])**2).sum(-1))
    gd=np.abs(psub[:,None]-psub[None,:])
    iu=np.triu_indices(len(ii),1)
    sp=stats.spearmanr(fd[iu], gd[iu]).statistic
    lab=AgglomerativeClustering(n_clusters=4).fit_predict(A)
    order=np.argsort(pos); sizes=np.bincount(lab, minlength=4)
    blk=np.empty(n, dtype=int); s=0
    for b,sz in enumerate(sizes):
        blk[order[s:s+sz]]=b; s+=sz
    ari=adjusted_rand_score(lab, blk)
    mb=np.digitize(mafv, [0.001,0.01,0.05])
    mi=0.0
    if len(np.unique(mb))>1 and len(np.unique(lab))>1:
        ct=pd.crosstab(lab, mb).values.astype(float); p=ct/ct.sum()
        px=p.sum(1, keepdims=True); py=p.sum(0, keepdims=True)
        with np.errstate(divide="ignore", invalid="ignore"):
            mi=float(np.nansum(p*np.log(p/(px*py))))
    rows.append(dict(domain=g, n=n, spearman_featdist_vs_genodist=float(sp),
                     ari_vs_position_blocks=float(ari), mi_cluster_vs_mafbin=mi))
pos_audit=pd.DataFrame(rows)
pos_audit.to_csv(f"{OUT}/prelim_position_audit.csv", index=False)
f4={"n_domains":int(len(pos_audit)),
    "spearman_median":float(pos_audit.spearman_featdist_vs_genodist.median()),
    "spearman_iqr":[float(pos_audit.spearman_featdist_vs_genodist.quantile(.25)), float(pos_audit.spearman_featdist_vs_genodist.quantile(.75))],
    "ari_position_median":float(pos_audit.ari_vs_position_blocks.median()),
    "ari_position_iqr":[float(pos_audit.ari_vs_position_blocks.quantile(.25)), float(pos_audit.ari_vs_position_blocks.quantile(.75))],
    "frac_domains_ari_ge_0.5":float((pos_audit.ari_vs_position_blocks>=0.5).mean()),
    "mi_maf_median":float(pos_audit.mi_cluster_vs_mafbin.median())}
verdict={"tendency_above_shuffle": bool(f1["hopkins_real_ci"][0] > f1["hopkins_shuffled_ci"][1]),
         "position_rediscovery_low": bool(f4["ari_position_median"] < 0.5)}
verdict["minimum_viability"]="PLAUSIBLE" if all(verdict.values()) else "NOT_SUPPORTED"
json.dump({"f1":f1,"f4_prelim":f4,"verdict":verdict,
           "vset_columns_untouched":vcols,"chroms":sorted(df.chr.unique().tolist())},
          open(f"{OUT}/prelim_verdict.json","w"), indent=1)
print(json.dumps({"f1":f1,"f4":f4,"verdict":verdict}, indent=1)[:2600])
