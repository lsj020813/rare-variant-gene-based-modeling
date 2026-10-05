
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import glob, json, csv
rows=[]; meta=None
for f in sorted(glob.glob(_config_path("${PROJECT_ROOT}/work/gate2/out/g2_chr*.json"))):
    d=json.load(open(f)); meta={k:d[k] for k in ("seed","nfold","nperm","maxpair","min_cocarrier","lambda_grid","N","primary_metric")}
    rows += d["genes"]
keys=sorted({k for r_ in rows for k in r_ if not isinstance(r_[k], list)})
with open(_config_path("${PROJECT_ROOT}/work/gate2/out/gate2_per_gene.csv"),"w",newline="") as fh:
    w=csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore"); w.writeheader()
    for r_ in rows: w.writerow({k:r_.get(k) for k in keys})
json.dump(meta, open(_config_path("${PROJECT_ROOT}/work/gate2/out/gate2_meta.json"),"w"))
print(json.dumps({"genes":len(rows),"ok":sum(1 for r_ in rows if r_.get("status")=="ok"),"files":len(glob.glob(_config_path("${PROJECT_ROOT}/work/gate2/out/g2_chr*.json"))),"meta":meta}))
