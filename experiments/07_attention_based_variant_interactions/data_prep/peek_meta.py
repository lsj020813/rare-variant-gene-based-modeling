
import sys, numpy as np, collections
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "data_prep"))
from prs_common import *
Ls=loci(); m=meta(Ls[0]); print("keys", list(m[0].keys()))
cnt=collections.Counter()
for L in Ls:
    for x in meta(L):
        f=float(x["maf_train"]); f=min(f,1-f)
        b="c>=5" if f>=0.05 else "l1-5" if f>=0.01 else "r0.1-1" if f>=0.001 else "v<0.1"
        cnt[b]+=1
print(dict(cnt), sum(cnt.values()))

import numpy as np
r2s=collections.defaultdict(list)
for L in Ls:
    for x in meta(L):
        f=float(x["maf_train"]); f=min(f,1-f)
        b="c>=5" if f>=0.05 else "l1-5" if f>=0.01 else "r0.1-1" if f>=0.001 else "v<0.1"
        r2s[b].append(float(x["r2"]))
print({k:(len(v), round(float(np.mean(v)),3)) for k,v in r2s.items()})
