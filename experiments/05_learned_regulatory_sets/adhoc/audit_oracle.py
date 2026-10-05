
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import json, glob, csv, math
A=_config_path("${PROJECT_ROOT}/work/fset/att")
rows=[]
for f in glob.glob(A+"/oracle/oracle_chr*.json"):
    d=json.load(open(f)); rows += d["per_domain"]
print("n_rows", len(rows), "keys", list(rows[0].keys())[:8])
bbj=set()
for r in csv.DictReader(open(A+"/att_bbj_tc_domains.csv")):
    if r['bbj_tc_domain']=='1': bbj.add(int(r['domain_idx']))
print("bbj n", len(bbj))
def binom_z(k,n,p): return (k-n*p)/math.sqrt(n*p*(1-p))
for s in ["ALL","COMMON","LOW","RARE","VRARE"]:
    ps=[r[s]["oof_p"] for r in rows if isinstance(r.get(s),dict) and r[s].get("oof_p") is not None]
    vals=sorted(set(round(p,4) for p in ps))
    k=sum(p<0.1 for p in ps); n=len(ps)
    ex=[r[s]["oof_z"]-r[s]["oof_z_null_mean"] for r in rows if isinstance(r.get(s),dict) and r[s].get("oof_z") is not None]
    ssx=[r[s]["single_z"]-r[s]["single_z_null_mean"] for r in rows if isinstance(r.get(s),dict) and r[s].get("single_z") is not None]
    b=[r for r in rows if r["domain_idx"] in bbj and isinstance(r.get(s),dict) and r[s].get("oof_p") is not None]
    bk=sum(r[s]["oof_p"]<0.1 for r in b)
    bex=sorted(r[s]["oof_z"]-r[s]["oof_z_null_mean"] for r in b)
    print(s, "n",n,"min_p",min(ps),"distinct_p",vals[:3],"k(p<.1)",k,"exp",round(n/11,1),"z",round(binom_z(k,n,1/11),2),
          "mean_ex",round(sum(ex)/len(ex),3),"mean_single_ex",round(sum(ssx)/len(ssx),3),
          "| BBJ n",len(b),"k",bk,"exp",round(len(b)/11,1),"median_ex",round(bex[len(bex)//2],3) if bex else None)
