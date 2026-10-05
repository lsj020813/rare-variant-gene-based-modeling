
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import json
d=json.load(open(_config_path("${PROJECT_ROOT}/work/phi_gate/peek/chr12_teamsize_diag.json")))
print("size", d["size_distribution"]); print("D1", d["D1_size_vs_target"])
R=d["per_annotation"]
def agg(k):
    v=[r[k]["mean"] for r in R if r[k]["mean"] is not None]
    import statistics
    return dict(n=len(v), mean=round(statistics.mean(v),4), median=round(statistics.median(v),4), pos=sum(x>0 for x in v))
for k in ("original","size_vs_annot","singleton_only","partial_size"): print(k, agg(k))
R2=sorted([r for r in R if r["original"]["genes"]>=100], key=lambda r:-r["original"]["mean"])
for r in R2[:8]+R2[-3:]:
    print(r["annotation"], *[(k, round(r[k]["mean"],4) if r[k]["mean"] is not None else None, round(r[k]["se"],4) if r[k]["se"] else None, r[k]["genes"]) for k in ("original","size_vs_annot","singleton_only","partial_size")])
