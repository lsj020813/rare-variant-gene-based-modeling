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


import sys, json, os, resource, time
sys.path.insert(0, _config_path("${PROJECT_ROOT}/work/phi_gate"))
resource.setrlimit(resource.RLIMIT_AS, (32 * 1024**3, 32 * 1024**3))
for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[k] = "4"
import numpy as np
from scipy.stats import rankdata
import final_phi_gate_v9 as G
G.np = np
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import Ridge
SMOKE = len(sys.argv) > 1 and sys.argv[1] == "smoke"
SEED = 20260929; P = 50 if SMOKE else 1000
def log(**k): k["t"] = time.strftime("%H:%M:%S"); print(json.dumps(k), flush=True)

def load(ch):
    b = G.read_built(ch, G.implementation())
    old = list(map(str, b["names"]))
    gl = np.asarray(b["gene_lengths"]); starts = np.r_[0, np.cumsum(gl)]
    ml = np.asarray(b["member_lengths"]); offs = np.r_[0, np.cumsum(ml)]; flat = np.asarray(b["member_flat"])
    res, _ = G.residual_target(b["absbeta"], b["maf"], b["lead"])
    tm = np.full(len(ml), np.nan)
    for t in range(len(ml)):
        v = res[flat[offs[t]:offs[t+1]]]; v = v[np.isfinite(v)]
        if len(v): tm[t] = v.mean()
    gene = np.repeat(np.arange(len(gl)), gl)
    st = np.select([ml == 1, ml == 2, ml <= 4], [0, 1, 2], 3)
    return dict(X=b["X"], names=old, y=tm, gene=gene, st=st, size=ml)

def align(d, names):
    X = np.zeros((len(d["X"]), len(names)))
    for j, n in enumerate(names):
        if n in d["names"]: X[:, j] = d["X"][:, d["names"].index(n)]
        else:
            assert n.startswith(("ccre_class=", "rep_class=")), n
    return X

def keep(d, X):
    ok = np.isfinite(d["y"]); return X[ok], d["y"][ok], d["gene"][ok], d["st"][ok], d["size"][ok]

def cell_demean(y, g, s):
    key = g * 4 + s; out = y.copy()
    u, inv = np.unique(key, return_inverse=True)
    m = np.bincount(inv, weights=y) / np.bincount(inv)
    return y - m[inv]

def zr(v):
    r = rankdata(v); r = r - r.mean(); n = np.linalg.norm(r); return None if n == 0 else r / n

def evaluate(pred, y, g, s, rng, singleton=False):
    obs, null, used = 0.0, np.zeros(P), 0
    for gg in np.unique(g):
        idx = np.where(g == gg)[0]
        if singleton: idx = idx[s[idx] == 0]
        if len(idx) < 3: continue
        yz, pz = zr(y[idx]), zr(pred[idx])
        if yz is None or pz is None: continue
        st = s[idx]; perm = np.tile(np.arange(len(idx)), (P, 1))
        for k in np.unique(st):
            pos = np.where(st == k)[0]
            if len(pos) > 1: perm[:, pos] = pos[np.argsort(rng.random((P, len(pos))), axis=1)]
        obs += yz.dot(pz); null += yz[perm].dot(pz); used += 1
    o = obs / used; nl = null / used
    return dict(genes=used, observed=float(o), null_mean=float(nl.mean()), null_q95=float(np.quantile(nl, .95)),
                null_sd=float(nl.std()), p=float((1 + (nl >= o).sum()) / (P + 1)))

def gb(): return HistGradientBoostingRegressor(max_iter=300, learning_rate=0.05, max_leaf_nodes=31, min_samples_leaf=200,
                                              l2_regularization=1.0, early_stopping=False, random_state=SEED)

d12 = load("12")
if SMOKE:
    names = sorted(d12["names"]); X12 = align(d12, names)
    Xa, ya, ga, sa, _ = keep(d12, X12)
    half = np.median(np.unique(ga)); tr = ga < half; te = ~tr
    Xtr, ytr, gtr, s_tr = Xa[tr], ya[tr], ga[tr], sa[tr]; Xte, yte, gte, s_te = Xa[te], ya[te], ga[te], sa[te]
else:
    d2 = load("2")
    names = sorted(set(d12["names"]) | set(d2["names"]))
    Xtr, ytr, gtr, s_tr, _ = keep(d12, align(d12, names)); Xte, yte, gte, s_te, _ = keep(d2, align(d2, names))
log(stage="data", train_teams=len(ytr), train_genes=len(np.unique(gtr)), eval_teams=len(yte), eval_genes=len(np.unique(gte)), features=len(names))
ttr = cell_demean(ytr, gtr, s_tr)
rng = np.random.default_rng(SEED)
ug = np.unique(gtr); fold = {g: i % 5 for i, g in enumerate(np.random.default_rng(SEED).permutation(ug))}
f = np.array([fold[g] for g in gtr])
med = np.nanmedian(Xtr, 0); med = np.where(np.isfinite(med), med, 0.0)
def std_fit(X):
    Z = np.where(np.isfinite(X), X, med); mu = Z.mean(0); sd = Z.std(0); sd[sd == 0] = 1; return mu, sd
def std_apply(X, mu, sd): return (np.where(np.isfinite(X), X, med) - mu) / sd
oof_gb = np.full(len(ytr), np.nan); alphas = [0.1, 1, 10, 100, 1000]; mse = {a: 0.0 for a in alphas}
for k in range(5):
    a_, b_ = f != k, f == k
    oof_gb[b_] = gb().fit(Xtr[a_], ttr[a_]).predict(Xtr[b_])
    mu, sd = std_fit(Xtr[a_])
    for a in alphas:
        r = Ridge(alpha=a).fit(std_apply(Xtr[a_], mu, sd), ttr[a_]); mse[a] += float(((r.predict(std_apply(Xtr[b_], mu, sd)) - ttr[b_]) ** 2).sum())
    log(stage="cv_fold", fold=k)
best_a = min(alphas, key=lambda a: mse[a])
oof = evaluate(oof_gb, ytr, gtr, s_tr, np.random.default_rng(SEED + 1))
log(stage="chr12_oof_gb", **oof, ridge_alpha=best_a)
m = gb().fit(Xtr, ttr); p_gb = m.predict(Xte)
mu, sd = std_fit(Xtr); r = Ridge(alpha=best_a).fit(std_apply(Xtr, mu, sd), ttr); p_r = r.predict(std_apply(Xte, mu, sd))
main = evaluate(p_gb, yte, gte, s_te, rng)
res = dict(scope=("SMOKE chr12 half/half" if SMOKE else "ML gate v1: train chr12, evaluate chr2 once"),
           permutations=P, features=len(names), train_teams=int(len(ytr)), eval_teams=int(len(yte)),
           primary_M_GB=main, pass_=bool(main["observed"] > main["null_q95"] and main["p"] < 0.05),
           secondary=dict(M_GB_singleton=evaluate(p_gb, yte, gte, s_te, rng, singleton=True),
                          M_R=evaluate(p_r, yte, gte, s_te, rng), ridge_alpha=best_a, chr12_oof_M_GB=oof))
res["verdict"] = "ML_GATE_PASS" if res["pass_"] else "ML_GATE_FAIL"
os.makedirs(_config_path("${PROJECT_ROOT}/work/phi_gate/peek"), exist_ok=True)
fn = _config_path("${PROJECT_ROOT}/work/phi_gate/peek/") + ("ml_gate_SMOKE.json" if SMOKE else "ml_gate_v1.json")
json.dump(res, open(fn + ".tmp", "w"), indent=1); os.replace(fn + ".tmp", fn); log(stage="done", verdict=res["verdict"] if not SMOKE else "smoke")
