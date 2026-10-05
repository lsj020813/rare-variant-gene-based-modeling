#!/usr/bin/env python
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import os, sys, json
import pandas as pd

BASE = _config_path("${PROJECT_ROOT}/work/fset/settest")
PREP = os.path.join(BASE, "prep")
GP = os.path.join(BASE, "groupfiles"); os.makedirs(GP, exist_ok=True)

rows = []
for N in range(1, 23):
    f = os.path.join(PREP, "chr%d.vars.tsv" % N)
    v = pd.read_csv(f, sep="\t", dtype={"domain": str})
    v["vid"] = "%d:" % N + v.pos.astype(str) + ":" + v.ref + ":" + v.alt
    out = os.path.join(GP, "chr%d.set.txt" % N)
    with open(out, "w") as fh:
        for dom, g in v.groupby("domain", sort=True):
            sets = {"all": g, "ccre": g[g.ccre == 1],
                    "major": g[g.is_major == 1], "minor": g[g.is_major == 0]}
            for nm, gg in sets.items():
                nv = len(gg)
                rows.append({"chr": N, "domain": dom, "set": nm, "n_var": nv,
                             "forced": int(g.forced.iloc[0]),
                             "written": int(nv > 0)})
                if nv == 0:
                    continue
                ids = " ".join(gg.vid.tolist())
                fh.write("%s__%s var %s\n" % (dom, nm, ids))
                fh.write("%s__%s anno %s\n" % (dom, nm, " ".join(["all"] * nv)))
    print("chr%d groups=%d" % (N, sum(1 for r in rows if r["chr"] == N and r["written"])), flush=True)

man = pd.DataFrame(rows)
man.to_csv(os.path.join(PREP, "group_manifest.csv"), index=False)
summ = {"n_groups_written": int(man.written.sum()),
        "n_groups_empty": int((1 - man.written).sum()),
        "empty_by_set": man[man.written == 0].groupby("set").size().to_dict(),
        "n_domains": int(man.domain.nunique())}
json.dump(summ, open(os.path.join(PREP, "group_summary.json"), "w"), indent=1)
print(json.dumps(summ))
