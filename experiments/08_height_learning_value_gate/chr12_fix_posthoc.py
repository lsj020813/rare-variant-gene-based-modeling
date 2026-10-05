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
resource.setrlimit(resource.RLIMIT_AS, (24 * 1024**3, 24 * 1024**3))
for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[k] = "2"
import numpy as np
from scipy.stats import rankdata
import final_phi_gate_v9 as G
G.np = np
P = 1000; SEED = 20260929
b = G.read_built("12", G.implementation())
old = list(map(str, b["names"])); names = sorted(old)
X = np.column_stack([b["X"][:, old.index(n)] for n in names])
gl = np.asarray(b["gene_lengths"]); starts = np.r_[0, np.cumsum(gl)]
ml = np.asarray(b["member_lengths"]); offs = np.r_[0, np.cumsum(ml)]; flat = np.asarray(b["member_flat"])
first = G.measure_one([(X[starts[g]:starts[g+1]], None) for g in range(len(gl))], names)
cand = [j for j, r in enumerate(first) if r["candidate"]]
res, _ = G.residual_target(b["absbeta"], b["maf"], b["lead"])
tm = np.full(len(ml), np.nan)
for t in range(len(ml)):
    v = res[flat[offs[t]:offs[t+1]]]; v = v[np.isfinite(v)]
    if len(v): tm[t] = v.mean()
stratum = np.select([ml == 1, ml == 2, ml <= 4], [0, 1, 2], 3)
rng = np.random.default_rng(SEED)
def zr(v):
    r = rankdata(v); r = r - r.mean(); n = np.linalg.norm(r); return None if n == 0 else r / n
def run(single):
    obs = np.zeros(len(cand)); null = np.zeros((P, len(cand))); used = np.zeros(len(cand))
    for g in range(len(gl)):
        idx = np.arange(starts[g], starts[g+1])
        ok = np.isfinite(tm[idx])
        if single: ok &= (ml[idx] == 1)
        idx = idx[ok]
        if len(idx) < 3: continue
        y = tm[idx]; st = stratum[idx]
        perm = np.tile(np.arange(len(idx)), (P, 1))
        for s in np.unique(st):
            pos = np.where(st == s)[0]
            if len(pos) > 1:
                perm[:, pos] = pos[np.argsort(rng.random((P, len(pos))), axis=1)]
        Xg = X[idx][:, cand]
        full = np.isfinite(Xg).all(0)
        cols = np.where(full)[0]
        if len(cols):
            yz = zr(y)
            if yz is not None:
                Z = np.column_stack([zr(Xg[:, c]) if zr(Xg[:, c]) is not None else np.full(len(idx), np.nan) for c in cols])
                good = np.isfinite(Z).all(0)
                if good.any():
                    cg = cols[good]; Zg = Z[:, good]
                    obs[cg] += yz.dot(Zg); null[:, cg] += yz[perm].dot(Zg); used[cg] += 1
        for c in np.where(~full)[0]:
            v = np.isfinite(Xg[:, c])
            if v.sum() < 3: continue
            yz = zr(y[v]); xz = zr(Xg[v, c])
            if yz is None or xz is None: continue
            sub = np.where(v)[0]; remap = -np.ones(len(idx), int); remap[sub] = np.arange(len(sub))
            pv = np.tile(np.arange(len(sub)), (P, 1)); stv = st[sub]
            for s in np.unique(stv):
                pos = np.where(stv == s)[0]
                if len(pos) > 1: pv[:, pos] = pos[np.argsort(rng.random((P, len(pos))), axis=1)]
            obs[c] += yz.dot(xz); null[:, c] += yz[pv].dot(xz); used[c] += 1
    out = []
    for k, j in enumerate(cand):
        if used[k] == 0: out.append(dict(annotation=names[j], genes=0, p=1.0)); continue
        o = obs[k] / used[k]; nl = null[:, k] / used[k]
        out.append(dict(annotation=names[j], genes=int(used[k]), observed=float(o), null_mean=float(nl.mean()),
                        null_q95=float(np.quantile(nl, .95)), null_sd=float(nl.std()), p=float((1 + (nl >= o).sum()) / (P + 1))))
    p = np.array([r["p"] for r in out]); m = len(p); o = np.argsort(p)
    q = np.empty(m); q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    for r, qq in zip(out, q): r["q"] = float(min(qq, 1)); r["pass"] = bool(r.get("observed") is not None and r["observed"] > r["null_q95"] and r["q"] < 0.10)
    return out
main = run(False); print("MAIN_DONE", sum(r["pass"] for r in main), flush=True)
sing = run(True); print("SINGLE_DONE", sum(r["pass"] for r in sing), flush=True)
os.makedirs(_config_path("${PROJECT_ROOT}/work/phi_gate/peek"), exist_ok=True)
tmp = _config_path("${PROJECT_ROOT}/work/phi_gate/peek/chr12_fix_posthoc.json.tmp")
json.dump(dict(scope="chr12 only, POST-HOC, team mean target + size-stratified within-gene permutation",
               candidates=len(cand), teams_with_target=int(np.isfinite(tm).sum()), main=main, singleton=sing), open(tmp, "w"), indent=1)
os.replace(tmp, tmp[:-4]); print("FIX_DONE")
