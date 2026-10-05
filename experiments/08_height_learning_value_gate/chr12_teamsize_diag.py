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
import numpy as np
from scipy.stats import rankdata
import final_phi_gate_v9 as G
G.np = np
b = G.read_built("12", G.implementation())
old = list(map(str, b["names"])); names = sorted(old)
x = np.column_stack([b["X"][:, old.index(n)] for n in names])
starts = np.r_[0, np.cumsum(b["gene_lengths"])]
grouped = [(x[starts[g]:starts[g+1]], np.full(b["gene_lengths"][g], np.nan)) for g in range(len(b["gene_lengths"]))]
offsets = np.r_[0, np.cumsum(b["member_lengths"])]
perteam = [b["member_flat"][offsets[t]:offsets[t+1]] for t in range(len(b["member_lengths"]))]
members = [perteam[starts[g]:starts[g+1]] for g in range(len(b["gene_lengths"]))]
size = np.asarray(b["member_lengths"], float)
sizes = [size[starts[g]:starts[g+1]] for g in range(len(b["gene_lengths"]))]
first = G.measure_one(grouped, names)
target, _ = G.residual_target(b["absbeta"], b["maf"], b["lead"])
G.attach_team_targets(grouped, members, target)
cand = [j for j, r in enumerate(first) if r["candidate"]]
def sp(a, c):
    a, c = rankdata(a), rankdata(c); a = a - a.mean(); c = c - c.mean()
    s = np.linalg.norm(a) * np.linalg.norm(c)
    return None if s == 0 else float(a.dot(c) / s)
def resid(v, z):
    v = rankdata(v); z = rankdata(z); Z = np.column_stack([np.ones(len(z)), z])
    return v - Z.dot(np.linalg.lstsq(Z, v, rcond=None)[0])
def summ(vals):
    v = np.array([t for t in vals if t is not None]);
    return dict(genes=int(len(v)), mean=float(v.mean()) if len(v) else None,
                se=float(v.std(ddof=1) / np.sqrt(len(v))) if len(v) > 1 else None, frac_pos=float((v > 0).mean()) if len(v) else None)
d1 = []
for (xx, y), s in zip(grouped, sizes):
    ok = np.isfinite(y)
    if ok.sum() >= 3 and np.ptp(s[ok]) > 0: d1.append(sp(s[ok], y[ok]))
size_dist = dict(frac_singleton=float((size == 1).mean()), median=float(np.median(size)), q90=float(np.quantile(size, .9)), max=float(size.max()))
rows = []
for j in cand:
    orig, d2, d3, d4 = [], [], [], []
    for (xx, y), s in zip(grouped, sizes):
        ok = np.isfinite(xx[:, j]) & np.isfinite(y)
        if ok.sum() >= 3:
            orig.append(sp(xx[ok, j], y[ok]))
            if np.ptp(s[ok]) > 0:
                d2.append(sp(s[ok], xx[ok, j]))
                a, c = resid(xx[ok, j], s[ok]), resid(y[ok], s[ok])
                n = np.linalg.norm(a) * np.linalg.norm(c); d4.append(None if n == 0 else float(a.dot(c) / n))
        one = ok & (s == 1)
        if one.sum() >= 3: d3.append(sp(xx[one, j], y[one]))
    rows.append(dict(annotation=names[j], original=summ(orig), size_vs_annot=summ(d2), singleton_only=summ(d3), partial_size=summ(d4)))
os.makedirs(_config_path("${PROJECT_ROOT}/work/phi_gate/peek"), exist_ok=True)
out = dict(scope="chr12 only, descriptive diagnostic, NOT the sealed verdict", size_distribution=size_dist,
           D1_size_vs_target=summ(d1), per_annotation=rows)
tmp = _config_path("${PROJECT_ROOT}/work/phi_gate/peek/chr12_teamsize_diag.json.tmp")
json.dump(out, open(tmp, "w"), indent=1); os.replace(tmp, tmp[:-4]); print("DIAG_DONE", len(rows))
