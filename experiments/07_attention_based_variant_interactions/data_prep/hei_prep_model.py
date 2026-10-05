#!/usr/bin/env python
import os as _os
if not _os.environ.get("PROJECT_ROOT", "").strip():
    raise ValueError("PROJECT_ROOT must be set and nonblank")
import numpy as np, csv, glob, json
W=(_os.environ["PROJECT_ROOT"] + '/work/prs'); P=W+"/private/hei"
sp=[l.rstrip("\n").split("\t") for l in open(P+"/split_hei.tsv")][1:]
ids=[a for a,_ in sp]; spl=np.array([b for _,b in sp]); tr=spl=="train"; va=spl=="valid"
ph={r["sample_id"]:[float(r[k]) for k in (["y"] + _os.environ["PHENO_COVARIATE_COLUMNS"].split(","))] for r in csv.DictReader(open((_os.environ["PROJECT_ROOT"] + '/work/ref/pheno_hei/hei.tsv')),delimiter="\t")}
A=np.array([ph[s] for s in ids]); y=A[:,0]; X=np.column_stack([np.ones(len(y)),A[:,1:]])
b=np.linalg.lstsq(X[tr],y[tr],rcond=None)[0]; r=y-X@b; r=r/r[tr].std()
sc={}
for f in sorted(glob.glob(W+"/score_hei/chr*.sscore")):
    with open(f) as g:
        h=g.readline().lstrip("#").split(); iid=h.index("IID"); col=h.index("SCORE1_SUM")
        for l in g:
            p=l.split(); sc[p[iid]]=sc.get(p[iid],0.0)+float(p[col])
assert len(glob.glob(W+"/score_hei/chr*.sscore"))==22
prs=np.array([sc.get(s,np.nan) for s in ids]); assert np.isfinite(prs).all()
prs=(prs-prs[tr].mean())/prs[tr].std()
np.savez(P+"/base_hei.npz", r=r.astype(np.float32), prs=prs.astype(np.float32), split=spl)
def r2(a,c): return float(1-((a-c)**2).sum()/((a-a.mean())**2).sum())
bp=np.polyfit(prs[tr],r[tr],1)
out=dict(n=len(ids),n_train=int(tr.sum()),cov_r2_train_raw=r2(y[tr],(X@b)[tr]),prs_r2_train=r2(r[tr],np.polyval(bp,prs[tr])),prs_r2_valid=r2(r[va],np.polyval(bp,prs[va])))
json.dump(out,open(W+"/out/hei/prep_model_hei.json","w"),indent=1); print(json.dumps(out))
