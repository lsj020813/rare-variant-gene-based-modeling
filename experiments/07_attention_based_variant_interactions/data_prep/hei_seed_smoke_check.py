#!/usr/bin/env python
import os as _os
if not _os.environ.get("PROJECT_ROOT", "").strip():
    raise ValueError("PROJECT_ROOT must be set and nonblank")
import numpy as np, json
W=(_os.environ["PROJECT_ROOT"] + '/work/prs')
def meta(f): return [l.split("\t")[0] for l in open(f).readlines()[1:]]
out={}
sk=meta(W+"/out/hei_seed_smoke/meta/meta_chr22.tsv"); SG=np.load(W+"/private/hei_seed_smoke/geno_chr22.npy")
ok=meta(W+"/out/hei/meta/meta_chr22.tsv"); OG=np.load(W+"/private/hei/geno_chr22.npy",mmap_mode="r")
oi={k:i for i,k in enumerate(ok)}; ov=[(i,oi[k]) for i,k in enumerate(sk) if k in oi]
out["dosage"]=dict(shape=list(SG.shape),n_keys=len(sk),overlap_old=len(ov),min=float(SG.min()),max=float(SG.max()),nan=int(np.isnan(SG.astype(np.float32)).sum()),
  max_abs_diff_vs_old=float(max([np.abs(SG[:,a].astype(np.float32)-OG[:,b].astype(np.float32)).max() for a,b in ov] or [-1])))
sk=[l.split("\t")[0] for l in open(W+"/out/hei_seed_smoke/gt/k_chr22.tsv").readlines()[1:]]; SH=np.load(W+"/private/hei_seed_smoke/gt/h_chr22.npy")
ok=[l.split("\t")[0] for l in open(W+"/out/hei/gt/k_chr22.tsv").readlines()[1:]]; OH=np.load(W+"/private/hei/gt/h_chr22.npy",mmap_mode="r")
oi={k:i for i,k in enumerate(ok)}; ov=[(i,oi[k]) for i,k in enumerate(sk) if k in oi]
dk={k:i for i,k in enumerate(meta(W+"/out/hei_seed_smoke/meta/meta_chr22.tsv"))}; both=[(i,dk[k]) for i,k in enumerate(sk) if k in dk]
out["gt"]=dict(shape=list(SH.shape),overlap_old=len(ov),n_mismatch_cells_vs_old=int(sum((SH[:,a,:]!=OH[:,b,:]).sum() for a,b in ov)),
  values=sorted(set(np.unique(SH).tolist())),strand_vs_dosage_overlap=len(both),
  corr_strandsum_dosage_min=float(min([np.corrcoef(SH[:,a,:].sum(1),SG[:,b].astype(np.float32))[0,1] for a,b in both] or [float("nan")])))
print(json.dumps(out,indent=1))
