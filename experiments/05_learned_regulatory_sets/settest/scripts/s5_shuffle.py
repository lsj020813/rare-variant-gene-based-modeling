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
import numpy as np, pandas as pd

BASE = _config_path("${PROJECT_ROOT}/work/fset/settest")
PREP = os.path.join(BASE, "prep"); GP = os.path.join(BASE, "groupfiles")
tag = sys.argv[1]
idx = [int(x) for x in sys.argv[2].split(",")]

rows = []
for N in range(1, 23):
    v = pd.read_csv(os.path.join(PREP, "chr%d.vars.tsv" % N), sep="\t", dtype={"domain": str})
    v["vid"] = "%d:" % N + v.pos.astype(str) + ":" + v.ref + ":" + v.alt
    out = os.path.join(GP, "chr%d.%s.txt" % (N, tag))
    with open(out, "w") as fh:
        for i in idx:
            rng = np.random.RandomState(20260923 + i)
            for dom, g in v.groupby("domain", sort=True):
                n = len(g); nmaj = int(g.is_major.sum())
                perm = rng.permutation(n)
                lab = np.zeros(n, dtype=int); lab[perm[:nmaj]] = 1
                ids = g.vid.values
                for nm, sel in (("major", lab == 1), ("minor", lab == 0)):
                    k = int(sel.sum())
                    rows.append({"chr": N, "shuffle": i, "domain": dom, "set": nm, "n_var": k})
                    if k == 0:
                        continue
                    gid = "%s__s%02d_%s" % (dom, i, nm)
                    fh.write("%s var %s\n" % (gid, " ".join(ids[sel].tolist())))
                    fh.write("%s anno %s\n" % (gid, " ".join(["all"] * k)))
    print("chr%d %s written" % (N, tag), flush=True)

pd.DataFrame(rows).to_csv(os.path.join(PREP, "shuffle_manifest_%s.csv" % tag), index=False)
print(json.dumps({"tag": tag, "shuffles": idx, "n_groups": len(rows)}))
