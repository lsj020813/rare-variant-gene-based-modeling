
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
B=_config_path("${PROJECT_ROOT}/work/phi_gate/")
m=json.load(open(B+"peek/ml_gate_v1.json"))
print("ML", m["verdict"], "train",m["train_teams"],"eval",m["eval_teams"])
for k,v in [("primary",m["primary_M_GB"]),("single",m["secondary"]["M_GB_singleton"]),("ridge",m["secondary"]["M_R"]),("oof12",m["secondary"]["chr12_oof_M_GB"])]:
    print(k, {a:(round(b,4) if isinstance(b,float) else b) for a,b in v.items()})
d=json.load(open(B+"out/pooled_2_12/decision.json"))
print("SEALED", d.get("decision"), "passed", len(d.get("passed_annotations",[])), "cand", d.get("candidates"))
f=json.load(open(B+"peek/chr2_fix_posthoc.json"))
for k in ("main","singleton"):
    R=[r for r in f[k] if r.get("observed") is not None]
    print("CHR2_FIX",k,"n",len(R),"pass",sum(r["pass"] for r in R),"p<.05",sum(r["p"]<.05 for r in R),"min_q",round(min(r["q"] for r in R),3),"mean_obs",round(statistics.mean(r["observed"] for r in R),4))
