
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import gzip, csv, collections, json, sys
import numpy as np
assign={}
algos=collections.Counter()
with gzip.open(_config_path("${PROJECT_ROOT}/work/fset/f1f2/ub_inub/consensus_clusters.csv.gz"),"rt") as f:
    rd=csv.DictReader(f)
    for r in rd:
        algos[(r["algo"],r["similarity"])]+=1
        if r["algo"]=="spectral" and r["similarity"]=="eucl":
            assign[(r["key"],r["domain"])]=r["consensus_cluster"]
bins=[(0,0.001,"<0.1%"),(0.001,0.01,"0.1-1%"),(0.01,0.05,"1-5%"),(0.05,0.2,"5-20%"),(0.2,0.51,">=20%")]
def b(m):
    for lo,hi,l in bins:
        if lo<=m<hi: return l
    return "NA"
cnt=collections.Counter(); tot=collections.Counter(); mafs=collections.defaultdict(list)
dom_rows=collections.defaultdict(list)
with open(_config_path("${PROJECT_ROOT}/work/fset/primary/primary_sample.tsv")) as f:
    rd=csv.DictReader(f, delimiter="\t")
    for r in rd:
        if r["in_ub"]!="1": continue
        k=(r["key"],r["domain"])
        if k in assign: dom_rows[r["domain"]].append((assign[k], float(r["maf"])))
n_dom=0
for d,rows in dom_rows.items():
    c=collections.Counter(a for a,_ in rows)
    if len(c)<2: continue
    n_dom+=1
    order=[a for a,_ in c.most_common()]
    for a,m in rows:
        lab="major" if a==order[0] else "minor"
        cnt[(lab,b(m))]+=1; tot[lab]+=1; mafs[lab].append(m)
out={"algos":{f"{a}|{s}":n for (a,s),n in algos.items()},"n_domains_with_2clusters":n_dom,"total_rows":{k:v for k,v in tot.items()}}
for lab in ("major","minor"):
    arr=np.array(mafs[lab])
    out[lab]={"n":len(arr),"median_maf":float(np.median(arr)),"q1":float(np.quantile(arr,.25)),"q3":float(np.quantile(arr,.75)),
              "bins":{l:cnt[(lab,l)] for _,_,l in bins},"bins_frac":{l:round(cnt[(lab,l)]/max(tot[lab],1),4) for _,_,l in bins}}
lower=0; nd=0
for d,rows in dom_rows.items():
    c=collections.Counter(a for a,_ in rows)
    if len(c)<2: continue
    order=[a for a,_ in c.most_common()]
    mj=np.median([m for a,m in rows if a==order[0]]); mn=np.median([m for a,m in rows if a!=order[0]])
    nd+=1; lower+= (mn<mj)
out["frac_domains_minor_cluster_lower_median_maf"]=round(lower/max(nd,1),4)
out["minor_cluster_size_frac_median"]=float(np.median([ (lambda c: c.most_common()[-1][1]/sum(c.values()))(collections.Counter(a for a,_ in rows)) for rows in dom_rows.values() if len(set(a for a,_ in rows))>=2]))
print(json.dumps(out,indent=1))
