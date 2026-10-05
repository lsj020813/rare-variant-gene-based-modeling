#!/usr/bin/env python
import os as _os
if not _os.environ.get("PROJECT_ROOT", "").strip():
    raise ValueError("PROJECT_ROOT must be set and nonblank")
import numpy as np, csv
W=(_os.environ["PROJECT_ROOT"] + '/work/prs'); P=W+"/private"
sp=[l.rstrip("\n").split("\t") for l in open(P+"/split.tsv")][1:]
ids=[a for a,_ in sp]; spl=np.array([b for _,b in sp]); tr=spl=="train"
ph={r["sample_id"]:[float(r[k]) for k in (["y"] + _os.environ["PHENO_COVARIATE_COLUMNS"].split(","))] for r in csv.DictReader(open((_os.environ["PROJECT_ROOT"] + '/work/ref/pheno_v3/tchl_v3.tsv')),delimiter="\t")}
A=np.array([ph[s] for s in ids]); y=A[:,0]; X=np.column_stack([np.ones(len(y)),A[:,1:]])
bb=np.linalg.lstsq(X[tr],y[tr],rcond=None)[0]; r=y-X@bb; r=r/r[tr].std()
np.savez(P+"/model/base_bench.npz", r=r.astype(np.float32), prs=np.zeros(len(r),dtype=np.float32), split=spl); print("ok")
