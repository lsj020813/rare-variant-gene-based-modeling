
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
    for r in csv.DictReader(f, delimiter="\t"):
        if r["in_ub"]!="1": continue
        k=(r["key"],r["domain"])
        if k in assign: dom[r["domain"]].append((assign[k], float(r["maf"])))
edges=np.linspace(-5,-0.3,48)
H={"major":np.zeros(len(edges)-1,int),"minor":np.zeros(len(edges)-1,int)}
permod=[]
for d,rows in dom.items():
    c=collections.Counter(a for a,_ in rows)
    if len(c)<2: continue
    top=c.most_common()[0][0]
    for lab,sel in (("major",lambda a:a==top),("minor",lambda a:a!=top)):
        m=np.log10(np.clip([x for a,x in rows if sel(a)],1e-5,0.5))
        H[lab]+=np.histogram(m,bins=edges)[0]
    permod.append((np.median(np.log10(np.clip([x for a,x in rows if a==top],1e-5,.5))), np.median(np.log10(np.clip([x for a,x in rows if a!=top],1e-5,.5)))))
print(json.dumps({"edges":edges.tolist(),"major":H["major"].tolist(),"minor":H["minor"].tolist(),"per_domain_median_log10maf":permod}))
