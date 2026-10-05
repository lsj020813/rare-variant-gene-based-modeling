import os as _os
if not _os.environ.get("PROJECT_ROOT", "").strip():
    raise ValueError("PROJECT_ROOT must be set and nonblank")

import csv, numpy as np
n=np.array([int(r["n_tokens"]) for r in csv.DictReader(open((_os.environ["PROJECT_ROOT"] + '/work/prs/out/hei/token_per_gene.tsv')),delimiter="\t")])
print("genes",len(n),"total",n.sum())
for cap in (64,128,256,512):
    m=np.minimum(n,cap); print(cap,"tokens",int(m.sum()),"attn_cost_rel_TCHL",round(float((m.astype(float)**2).sum()/(42*1024**2)),2))
import collections
L=collections.Counter(r["locus"] for r in csv.DictReader(open((_os.environ["PROJECT_ROOT"] + '/work/prs/out/hei/token_per_gene.tsv')),delimiter="\t"))
v=np.array(list(L.values())); print("genes_per_locus q",[float(np.percentile(v,q)) for q in (0,25,50,75,100)])
