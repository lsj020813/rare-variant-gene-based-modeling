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


import os, json
import pandas as pd, numpy as np

BASE = _config_path("${PROJECT_ROOT}/work/fset/settest")
PREP = os.path.join(BASE, "prep"); os.makedirs(PREP, exist_ok=True)
PS = _config_path("${PROJECT_ROOT}/work/fset/primary/primary_sample.tsv")
CC = _config_path("${PROJECT_ROOT}/work/fset/f1f2/ub_inub/consensus_clusters.csv.gz")

ps = pd.read_csv(PS, sep="\t", dtype={"key": str, "chr": str, "domain": str})
ps = ps[ps.in_ub == 1].copy()
assert len(ps) == 364430, len(ps)

cc = pd.read_csv(CC, dtype={"key": str, "domain": str})
cc = cc[(cc.algo == "spectral") & (cc.similarity == "eucl")][["key", "domain", "consensus_cluster"]]
assert len(cc) == 364430, len(cc)

m = ps.merge(cc, on=["key", "domain"], how="left", validate="one_to_one")
assert m.consensus_cluster.notna().all(), "unjoined rows"
m["consensus_cluster"] = m.consensus_cluster.astype(int)

sz = m.groupby(["domain", "consensus_cluster"]).size().rename("n").reset_index()
sz = sz.sort_values(["domain", "n", "consensus_cluster"], ascending=[True, False, True])
maj = sz.groupby("domain", sort=False).head(1)[["domain", "consensus_cluster"]].copy()
maj["is_major"] = 1
m = m.merge(maj, on=["domain", "consensus_cluster"], how="left")
m["is_major"] = m.is_major.fillna(0).astype(int)
kpd = sz.groupby("domain").size()

kp = m.key.str.split(":", expand=True)
m["vchr"] = kp[0].str.replace("^chr", "", regex=True)
m["pos"] = kp[1].astype(np.int64); m["ref"] = kp[2]; m["alt"] = kp[3]
assert (m.vchr == m["chr"].astype(str)).all(), "chr mismatch between key and chr column"
m["ccre"] = (m.has_ccre > 0).astype(int)

dom_forced = m.groupby("domain").forced.nunique()
assert (dom_forced == 1).all(), "forced flag not constant within domain"

stats = {
    "n_var": int(len(m)), "n_domain": int(m.domain.nunique()),
    "n_ccre_var": int(m.ccre.sum()), "n_major_var": int(m.is_major.sum()),
    "n_minor_var": int((1 - m.is_major).sum()),
    "n_domains_K2": int((kpd == 2).sum()), "n_domains_Kgt2": int((kpd > 2).sum()),
    "n_domains_forced": int(m.groupby("domain").forced.first().sum()),
    "n_domains_random": int((m.groupby("domain").forced.first() == 0).sum()),
    "n_chr": int(m.vchr.nunique()),
}

for ch, g in m.groupby("vchr", sort=False):
    g = g.sort_values("pos")
    g[["pos", "ref", "alt", "domain", "ccre", "is_major", "maf", "r2", "typed", "forced"]].to_csv(
        os.path.join(PREP, "chr%s.vars.tsv" % ch), sep="\t", index=False)
    with open(os.path.join(PREP, "chr%s.pos.bed" % ch), "w") as f:
        for p in np.unique(g.pos.values):
            f.write("%s\t%d\t%d\n" % (ch, p - 1, p))

json.dump(stats, open(os.path.join(PREP, "prep_stats.json"), "w"), indent=1)
print(json.dumps(stats))
