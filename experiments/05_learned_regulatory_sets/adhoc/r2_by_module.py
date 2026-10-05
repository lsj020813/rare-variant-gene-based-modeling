
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import gzip, csv, collections, json
import numpy as np
assign={}
with gzip.open(_config_path("${PROJECT_ROOT}/work/fset/f1f2/ub_inub/consensus_clusters.csv.gz"),"rt") as f:
    for r in csv.DictReader(f):
        if r["algo"]=="spectral" and r["similarity"]=="eucl": assign[(r["key"],r["domain"])]=r["consensus_cluster"]
dom=collections.defaultdict(list)
with open(_config_path("${PROJECT_ROOT}/work/fset/primary/primary_sample.tsv")) as f:
    rd=csv.DictReader(f, delimiter="\t"); cols=rd.fieldnames
    for r in rd:
        if r["in_ub"]!="1": continue
        k=(r["key"],r["domain"])
        if k in assign: dom[r["domain"]].append((assign[k], float(r["r2"]) if r["r2"] not in ("","NA") else np.nan, r["typed"]))
def r2bin(x):
    if np.isnan(x): return "NA"
    return "<0.3" if x<0.3 else "0.3-0.8" if x<0.8 else "0.8-0.9" if x<0.9 else ">=0.9"
cnt=collections.Counter(); tot=collections.Counter(); typed=collections.Counter(); r2s=collections.defaultdict(list)
per=[]
for d,rows in dom.items():
    c=collections.Counter(a for a,_,_ in rows)
    if len(c)<2: continue
    top=c.most_common()[0][0]
    for a,r2,t in rows:
        lab="major" if a==top else "minor"
        tot[lab]+=1; cnt[(lab,r2bin(r2))]+=1; typed[lab]+= (t in ("1","True","TYPED"))
        if not np.isnan(r2): r2s[lab].append(r2)
    mj=np.nanmedian([r2 for a,r2,_ in rows if a==top]); mn=np.nanmedian([r2 for a,r2,_ in rows if a!=top]); per.append((mj,mn))
per=np.array(per)
out={"typed_values_seen":list(collections.Counter(t for rows in dom.values() for _,_,t in rows).items())[:5]}
for lab in ("major","minor"):
    a=np.array(r2s[lab])
    out[lab]={"n":tot[lab],"r2_median":float(np.median(a)),"r2_q1":float(np.quantile(a,.25)),"r2_q3":float(np.quantile(a,.75)),
              "r2_bins_frac":{b:round(cnt[(lab,b)]/tot[lab],4) for b in ("<0.3","0.3-0.8","0.8-0.9",">=0.9","NA")},
              "typed_frac":round(typed[lab]/tot[lab],4)}
out["frac_domains_minor_lower_median_r2"]=float(np.mean(per[:,1]<per[:,0])); out["n_domains"]=len(per)
out["median_domain_r2_diff_minor_minus_major"]=float(np.median(per[:,1]-per[:,0]))
print(json.dumps(out,indent=1))
