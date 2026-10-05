
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
d=json.load(open(_config_path("${PROJECT_ROOT}/work/phi_gate/peek/chr12_fix_posthoc.json")))
print("cand",d["candidates"],"teams_w_target",d["teams_with_target"])
for k in ("main","singleton"):
    R=[r for r in d[k] if r.get("observed") is not None]
    v=[r["observed"] for r in R]; z=[(r["observed"]-r["null_mean"])/r["null_sd"] for r in R if r["null_sd"]>0]
    print(k,"n",len(R),"pass",sum(r["pass"] for r in R),"mean_obs",round(statistics.mean(v),4),"pos",sum(x>0 for x in v),"min_q",round(min(r["q"] for r in R),3),"n_p<.05",sum(r["p"]<.05 for r in R),"mean_z",round(statistics.mean(z),2))
    for r in sorted(R,key=lambda r:r["p"])[:5]: print("  ",r["annotation"],round(r["observed"],4),"null95",round(r["null_q95"],4),"p",round(r["p"],3),"q",round(r["q"],3),r["genes"])
