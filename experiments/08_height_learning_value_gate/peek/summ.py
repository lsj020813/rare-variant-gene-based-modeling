
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import json, statistics
d=json.load(open(_config_path("${PROJECT_ROOT}/work/phi_gate/peek/chr12_peek.json")))
m1=d["measurement_1"]
vr=sorted(r["variance_ratio"] for r in m1 if r["variance_ratio"] is not None)
print("m1 ratio q0,25,50,75,100", [round(vr[int(q*(len(vr)-1))],3) for q in (0,.25,.5,.75,1)])
m2=[r for r in d["measurement_2_descriptive"] if r["observed"] is not None]
m2.sort(key=lambda r:-r["observed"])
v=[r["observed"] for r in m2]
print("m2 n",len(v),"median",round(statistics.median(v),4),"mean",round(statistics.mean(v),4),"n>0",sum(x>0 for x in v),"n>0.02",sum(x>0.02 for x in v))
for r in m2[:10]: print("TOP",r["annotation"],round(r["observed"],4),r["eligible_genes"])
for r in m2[-5:]: print("BOT",r["annotation"],round(r["observed"],4),r["eligible_genes"])
print({k:v for k,v in d["regression"].items() if not isinstance(v,(list,dict))})
