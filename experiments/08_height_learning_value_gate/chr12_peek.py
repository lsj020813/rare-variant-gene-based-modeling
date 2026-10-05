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


import sys, json, os, resource
sys.path.insert(0, _config_path("${PROJECT_ROOT}/work/phi_gate"))
resource.setrlimit(resource.RLIMIT_AS, (16 * 1024**3, 16 * 1024**3))
for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[k] = "1"
import numpy
import final_phi_gate_v9 as G
G.np = numpy; np = numpy
b = G.read_built("12", G.implementation())
names = sorted(map(str, b["names"]))
old = list(map(str, b["names"]))
x = np.column_stack([b["X"][:, old.index(n)] for n in names])
starts = np.r_[0, np.cumsum(b["gene_lengths"])]
grouped = [(x[starts[g]:starts[g+1]], np.full(b["gene_lengths"][g], np.nan)) for g in range(len(b["gene_lengths"]))]
offsets = np.r_[0, np.cumsum(b["member_lengths"])]
perteam = [b["member_flat"][offsets[t]:offsets[t+1]] for t in range(len(b["member_lengths"]))]
members = [perteam[starts[g]:starts[g+1]] for g in range(len(b["gene_lengths"]))]
first = G.measure_one(grouped, names)
target, regression = G.residual_target(b["absbeta"], b["maf"], b["lead"])
G.attach_team_targets(grouped, members, target)
desc = G.descriptive_correlations(grouped, names, first)
os.makedirs(_config_path("${PROJECT_ROOT}/work/phi_gate/peek"), exist_ok=True)
res = dict(scope="chr12 only, descriptive, NOT the sealed verdict", genes=len(grouped),
           team_instances=int(sum(b["gene_lengths"])), regression=regression,
           measurement_1=first, measurement_2_descriptive=desc)
tmp = _config_path("${PROJECT_ROOT}/work/phi_gate/peek/chr12_peek.json.tmp")
json.dump(res, open(tmp, "w"), indent=1, default=float); os.replace(tmp, tmp[:-4])
print("PEEK_DONE", len(first), sum(1 for r in first if r.get("candidate")))
