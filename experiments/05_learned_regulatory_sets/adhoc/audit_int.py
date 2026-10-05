
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import json, glob, csv
from math import comb
A=_config_path("${PROJECT_ROOT}/work/fset/att")
bbj=set(int(r['domain_idx']) for r in csv.DictReader(open(A+"/att_bbj_tc_domains.csv")) if r['bbj_tc_domain']=='1')
def bt(k,n,p=1/11): return sum(comb(n,i)*p**i*(1-p)**(n-i) for i in range(k,n+1))
rr=[]
for f in glob.glob(A+"/interact/int_chr*.json"): rr+=json.load(open(f))["per_domain"]
for s in ["MAIN","MAIN_INT"]:
    a=[r[s] for r in rr if r.get(s)]; b=[r[s] for r in rr if r.get(s) and r["domain_idx"] in bbj]
    ka=sum(x["p"]<0.1 for x in a); kb=sum(x["p"]<0.1 for x in b)
    eb=sorted(x["excess"] for x in b)
    print(s,"all",len(a),ka,"exp",round(len(a)/11,1),"tail",round(bt(ka,len(a)),4),"| BBJ",len(b),kb,"tail",round(bt(kb,len(b)),6),"median_ex",round(eb[len(eb)//2],3))
d=sorted(r["MAIN_INT"]["excess"]-r["MAIN"]["excess"] for r in rr if r.get("MAIN") and r.get("MAIN_INT") and r["domain_idx"] in bbj)
print("BBJ INT-MAIN median",round(d[len(d)//2],3),"frac>0",round(sum(x>0 for x in d)/len(d),2), len(d))
