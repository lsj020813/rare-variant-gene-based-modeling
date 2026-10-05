#!/usr/bin/env python
import os as _os
if not _os.environ.get("PROJECT_ROOT", "").strip():
    raise ValueError("PROJECT_ROOT must be set and nonblank")
import csv, numpy as np, json
from scipy.stats import rankdata, norm
R=list(csv.DictReader(open((_os.environ["PROJECT_ROOT"] + '/work/ref/pheno_hei/hei.tsv')),delimiter="\t")); h=list(R[0].keys())
y=np.array([float(r["y"]) for r in R]); z=norm.ppf((rankdata(y)-0.375)/(len(y)+0.25))
with open((_os.environ["PROJECT_ROOT"] + '/work/ref/pheno_hei/hei_rint.tsv'),"w") as g:
    g.write("\t".join(h)+"\n")
    for r,v in zip(R,z): r=dict(r); r["y"]=f"{v:.6f}"; g.write("\t".join(r[k] for k in h)+"\n")
print(json.dumps(dict(n=len(y),cols=h)))
