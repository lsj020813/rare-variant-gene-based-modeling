
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import glob, json
import numpy as np, pandas as pd
ROOT=_config_path("${PROJECT_ROOT}/work"); OUT=f"{ROOT}/fset/prelim"; FO=f"{ROOT}/fset/out"
df=pd.concat([pd.read_csv(f, sep="\t", low_memory=False) for f in sorted(glob.glob(f"{FO}/feat_chr*.tsv.gz"))], ignore_index=True)
cat=["ccre_class","rep_class"]; 
tfc=[c for c in df.columns if c.startswith("tf_") and c!="tf_n"]
disc=df[cat].fillna("none").astype(str).agg("|".join, axis=1) + "|" + df["cpg"].fillna(0).astype(int).astype(str) + "|" + df[tfc].fillna(0).astype(int).astype(str).agg("".join, axis=1)
u=disc.value_counts()
num=["ccre_dist","rep_dist","cpg_dist","map_k36","tf_n"]
full=disc + "|" + df[num].round(3).astype(str).agg("|".join, axis=1)
uf=full.value_counts()
res={"n_variants":int(len(df)),
 "unique_discrete_vectors":int(u.size), "top_discrete_share":float(u.iloc[0]/len(df)),
 "share_in_discrete_groups_ge10":float(u[u>=10].sum()/len(df)),
 "share_singleton_discrete":float((u==1).sum()/len(df)),
 "unique_full_vectors":int(uf.size), "share_singleton_full":float((uf==1).sum()/len(df)),
 "ccre_none_frac":float((df.ccre_class.fillna("none")=="none").mean()),
 "tf_n_median":float(df.tf_n.median()), "tf_n_zero_frac":float((df.tf_n==0).mean()),
 "top5_discrete_groups":[[k[:60], int(v)] for k,v in u.head(5).items()]}
json.dump(res, open(f"{OUT}/prelim_duplicate_audit.json","w"), indent=1)
print(json.dumps(res, indent=1)[:1600])
