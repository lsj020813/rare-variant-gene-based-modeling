import os as _os
if not _os.environ.get("PROJECT_ROOT", "").strip():
    raise ValueError("PROJECT_ROOT must be set and nonblank")

import numpy as np, glob, os, gzip, csv, json
W=(_os.environ["PROJECT_ROOT"] + '/work/prs'); P=W+"/private"
def loci():
    return [l.split()[3] for l in open(W+"/out/loci.bed") if os.path.exists(f"{P}/loci/{l.split()[3]}.npz")]
def base(tag="bbj"):
    d=np.load(f"{P}/model/base_{tag}.npz", allow_pickle=True); return d["r"].astype(np.float64), d["prs"].astype(np.float64), d["split"]
def meta(L):
    rows=list(csv.DictReader(open(f"{W}/out/loci_meta/{L}.tsv"),delimiter="\t")); return rows
def dosage(L, tr):
    D=np.load(f"{P}/loci/{L}.npz")["D"].astype(np.float32)
    mu=D[tr].mean(0); sd=D[tr].std(0); sd[sd<1e-6]=1.0
    return (D-mu)/sd
def feats():
    F={}; cols=None
    for f in sorted(glob.glob(W+"/out/feat/feat_chr*.tsv.gz")):
        with gzip.open(f,"rt") as g:
            h=g.readline().rstrip("\n").split("\t")
            if cols is None:
                cols=[c for c in h if c not in ("key","chr","pos38","domain") and "re2g" not in c.lower() and not c.endswith("_target")]
            ix=[h.index(c) for c in cols]
            for l in g:
                p=l.rstrip("\n").split("\t"); F[p[0]]=[p[i] for i in ix]
    return F, cols
def encode_feats(F, cols, keys):
    num, cat = [], []
    for j,c in enumerate(cols):
        vals=[F.get(k,[""]*len(cols))[j] for k in keys[:5000]]
        try: [float(v) for v in vals if v!=""]; num.append(j)
        except ValueError: cat.append(j)
    out=[]; names=[]
    for j in num:
        v=np.array([float(F[k][j]) if (k in F and F[k][j]!="") else np.nan for k in keys])
        miss=np.isnan(v); v[miss]=0.0; out.append(v); names.append(cols[j])
        if miss.any() and not miss.all(): out.append(miss.astype(float)); names.append(cols[j]+"_miss")
    for j in cat:
        lv=sorted(set(F[k][j] for k in keys if k in F))
        for x in lv[:30]:
            out.append(np.array([1.0 if (k in F and F[k][j]==x) else 0.0 for k in keys])); names.append(f"{cols[j]}={x}")
    return np.column_stack(out) if out else np.zeros((len(keys),0)), names
def r2(y,p): return float(1-((y-p)**2).sum()/((y-y.mean())**2).sum())
