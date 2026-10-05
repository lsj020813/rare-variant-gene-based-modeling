
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
from math import comb
A=_config_path("${PROJECT_ROOT}/work/fset/att")
bbj=set(int(r['domain_idx']) for r in csv.DictReader(open(A+"/att_bbj_tc_domains.csv")) if r['bbj_tc_domain']=='1')
def btail(k,n,p=1/11): return sum(comb(n,i)*p**i*(1-p)**(n-i) for i in range(k,n+1))
rows=[]
for f in glob.glob(A+"/oracle/oracle_chr*.json"): rows+=json.load(open(f))["per_domain"]
for s in ["ALL","COMMON"]:
    b=[r[s] for r in rows if r["domain_idx"] in bbj and isinstance(r.get(s),dict) and r[s].get("oof_z") is not None]
    d=sorted((x["oof_z"]-x["oof_z_null_mean"])-(x["single_z"]-x["single_z_null_mean"]) for x in b)
    k=sum(x["oof_p"]<0.1 for x in b)
    print("ORACLE",s,"BBJ n",len(b),"k",k,"binom_tail",round(btail(k,len(b)),6),"median(oofEx-singleEx)",round(d[len(d)//2],3),"frac>0",round(sum(x>0 for x in d)/len(d),2))
for sub in ["oracle_hq","interact"]:
    fs=glob.glob(f"{A}/{sub}/*chr*.json")
    rr=[]
    for f in fs:
        try:
            d=json.load(open(f)); rr+=d.get("per_domain",[])
        except Exception as e: pass
    if not rr: print(sub,"no per_domain", fs[:2]); continue
    keys=[k for k in rr[0] if isinstance(rr[0][k],dict)]
    for s in keys:
        b=[r[s] for r in rr if r["domain_idx"] in bbj and r[s].get("oof_p") is not None]
        a=[r[s] for r in rr if r[s].get("oof_p") is not None]
        if not a: continue
        kb=sum(x["oof_p"]<0.1 for x in b); ka=sum(x["oof_p"]<0.1 for x in a)
        print(sub,s,"all n",len(a),"k",ka,"exp",round(len(a)/11,1),"| BBJ n",len(b),"k",kb,"tail",round(btail(kb,len(b)),6) if b else None)
