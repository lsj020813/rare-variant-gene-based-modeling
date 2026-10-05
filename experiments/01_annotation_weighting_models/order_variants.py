#!/usr/bin/env python3
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
O=_config_path("${PROJECT_ROOT}/work/run_trackB")
rows=[]
with open(f"{O}/variants.tsv") as fh:
    hdr=fh.readline()
    for line in fh: rows.append(line)
def key(l): t=l.split("\t"); return (-int(t[7]), int(t[2]))
rows.sort(key=key)
nt=sum(1 for l in rows if l.split("\t")[7]=="1")
with open(f"{O}/variants_ordered.tsv","w") as fo:
    fo.write(hdr); fo.writelines(rows)
json.dump({"n_total":len(rows),"n_target_first":nt,"n_ctrl_only":len(rows)-nt,"order":"is_target desc, pos asc"},open(f"{O}/variants_ordered_manifest.json","w"),indent=1)
print(json.dumps({"n_total":len(rows),"n_target_first":nt,"n_ctrl_only":len(rows)-nt}))
