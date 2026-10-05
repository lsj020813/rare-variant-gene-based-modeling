import os as _os
if not _os.environ.get("PROJECT_ROOT", "").strip():
    raise ValueError("PROJECT_ROOT must be set and nonblank")

import numpy as np, glob, csv
W=(_os.environ["PROJECT_ROOT"] + '/work/prs')
d=np.load(W+"/private/model/base_bbj.npz",allow_pickle=True); r=d["r"]; sp=d["split"]; tr=sp=="train"
ids=[l.split("\t")[0] for l in open(W+"/private/split.tsv")][1:]
pos={s:i for i,s in enumerate(ids)}
for N in [1,2,16,19,22]:
    v=np.full(len(ids),np.nan)
    with open(f"{W}/score_bbj/chr{N}.sscore") as g:
        g.readline()
        for l in g:
            p=l.split(); v[pos[p[0]]]=float(p[4])
    print(N, "corr_train", round(float(np.corrcoef(v[tr],r[tr])[0,1]),4), "sd", float(np.nanstd(v)))
import os
best=[]
for f in sorted(glob.glob(W+"/private/loci/L*.npz"))[:42]:
    D=np.load(f)["D"].astype(np.float32)[tr]
    Dc=D-D.mean(0); rc=r[tr]-r[tr].mean()
    cc=(Dc*rc[:,None]).sum(0)/np.sqrt((Dc**2).sum(0)*(rc**2).sum()+1e-12)
    best.append(float(np.nanmax(np.abs(cc))))
best=np.array(best); print("locus max|corr| median",round(float(np.median(best)),4),"max",round(float(best.max()),4),"n>0.03",int((best>0.03).sum()))
