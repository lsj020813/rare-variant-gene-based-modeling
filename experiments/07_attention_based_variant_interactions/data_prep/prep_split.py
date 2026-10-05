#!/usr/bin/env python
import os as _os
if not _os.environ["SPARSE_GRM_PATH"].strip():
    raise ValueError("Set nonblank SPARSE_GRM_PATH")
if not _os.environ.get("PROJECT_ROOT", "").strip():
    raise ValueError("PROJECT_ROOT must be set and nonblank")
import numpy as np, csv, json, collections, os
SEED=20260926
G=_os.environ["SPARSE_GRM_PATH"]
ids=[l.strip() for l in open(G+".sampleIDs.txt")]
ph={}
for r in csv.DictReader(open((_os.environ["PROJECT_ROOT"] + '/work/ref/pheno_v3/tchl_v3.tsv')), delimiter="\t"):
    try:
        v=[float(r[k]) for k in (["y"] + _os.environ["PHENO_COVARIATE_COLUMNS"].split(","))]
        if all(np.isfinite(v)): ph[r["sample_id"]]=1
    except: pass
def norm(s): return s.split("_")[0] if s not in ph and s.split("_")[0] in ph else s
pids=[norm(s) for s in ids]
match=sum(p in ph for p in pids)
edges=collections.defaultdict(set)
with open(G) as f:
    f.readline(); f.readline()
    for line in f:
        i,j,v=line.split(); i=int(i)-1; j=int(j)-1
        if i!=j:
            a,b=pids[i],pids[j]
            if a in ph and b in ph: edges[a].add(b); edges[b].add(a)
rng=np.random.default_rng(SEED)
removed=set(); deg={k:len(v) for k,v in edges.items()}
tieb={k:rng.random() for k in edges}
while True:
    live=[k for k in edges if k not in removed and any(n not in removed for n in edges[k])]
    if not live: break
    k=max(live, key=lambda x:(sum(n not in removed for n in edges[x]), tieb[x]))
    removed.add(k)
keep=sorted(p for p in ph if p not in removed)
rng2=np.random.default_rng(SEED); perm=rng2.permutation(len(keep))
n=len(keep); a=int(0.6*n); b=int(0.8*n)
lab=np.empty(n,dtype=object); lab[perm[:a]]="train"; lab[perm[a:b]]="valid"; lab[perm[b:]]="test"
os.makedirs((_os.environ["PROJECT_ROOT"] + '/work/prs/private'),exist_ok=True)
tmp=(_os.environ["PROJECT_ROOT"] + '/work/prs/private/split.tsv.tmp')
with open(tmp,"w") as f:
    f.write("sample_id\tsplit\n")
    for s,l in zip(keep,lab): f.write(f"{s}\t{l}\n")
os.replace(tmp,(_os.environ["PROJECT_ROOT"] + '/work/prs/private/split.tsv'))
out=dict(grm_ids=len(ids), pheno_complete=len(ph), grm_pheno_match=match, n_pairs_in_pheno=sum(len(v) for v in edges.values())//2,
         removed=len(removed), kept=n, train=int((lab=="train").sum()), valid=int((lab=="valid").sum()), test=int((lab=="test").sum()), seed=SEED)
json.dump(out, open((_os.environ["PROJECT_ROOT"] + '/work/prs/out/split_summary.json'),"w"), indent=1); print(json.dumps(out))
