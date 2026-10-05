#!/usr/bin/env python
import os as _os
if not _os.environ.get("PROJECT_ROOT", "").strip():
    raise ValueError("PROJECT_ROOT must be set and nonblank")
import numpy as np, csv, gzip, glob, os, json, subprocess, sys
W=(_os.environ["PROJECT_ROOT"] + '/work/prs'); P=W+"/private"; os.makedirs(P+"/model",exist_ok=True)
TAG=sys.argv[1] if len(sys.argv)>1 else "bbj"
sp=[l.rstrip("\n").split("\t") for l in open(P+"/split.tsv")][1:]
ids=[a for a,_ in sp]; spl=np.array([b for _,b in sp]); tr=spl=="train"
ph={}
for r in csv.DictReader(open((_os.environ["PROJECT_ROOT"] + '/work/ref/pheno_v3/tchl_v3.tsv')),delimiter="\t"):
    ph[r["sample_id"]]=[float(r[k]) for k in (["y"] + _os.environ["PHENO_COVARIATE_COLUMNS"].split(","))]
A=np.array([ph[s] for s in ids]); y=A[:,0]; X=np.column_stack([np.ones(len(y)),A[:,1:]])
b=np.linalg.lstsq(X[tr],y[tr],rcond=None)[0]; r=y-X@b; r=r/r[tr].std()
sc={}
for f in sorted(glob.glob(f"{W}/score_{TAG}/chr*.sscore")):
    with open(f) as g:
        h=g.readline().lstrip("#").split(); iid=h.index("IID"); col=h.index("SCORE1_SUM")
        for l in g:
            p=l.split(); sc[p[iid]]=sc.get(p[iid],0.0)+float(p[col])
prs=np.array([sc.get(s,np.nan) for s in ids]); assert np.isfinite(prs).all(), "missing prs"
prs=(prs-prs[tr].mean())/prs[tr].std()
np.savez(P+f"/model/base_{TAG}.npz", r=r.astype(np.float32), prs=prs.astype(np.float32), split=spl)
bp=np.polyfit(prs[tr], r[tr], 1)
def r2(a,b): return 1-((a-b)**2).sum()/((a-a.mean())**2).sum()
va=spl=="valid"
out=dict(tag=TAG, n=len(ids), n_train=int(tr.sum()), prs_r2_train=float(r2(r[tr],np.polyval(bp,prs[tr]))), prs_r2_valid=float(r2(r[va],np.polyval(bp,prs[va]))))
json.dump(out, open(W+f"/out/prep_model_{TAG}.json","w"), indent=1); print(json.dumps(out))
