
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import numpy as np, json, sys, time, os
NPERM = int(os.environ.get("NPERM", "50"))
sys.argv = ["x", "--fold", "0", "--arm", "spline", "--gradcheck-only"]
src = open("l1_train.py").read(); cut = src.index("if A.gcv_probe:")
g = {}; exec(compile(src[:cut], "l1_loader", "exec"), g)
TR, TRAITS, test_gi, burden, zrow, keys = g["TR"], g["TRAITS"], g["test_gi"], g["burden"], g["zrow"], g["keys"]
log = g["log"]
phi = np.load(_config_path("${PROJECT_ROOT}/work/ref/annot/l1/spline_fold0.v3.phi.npy")).astype(np.float32)
assert len(phi) == len(keys), (len(phi), len(keys))
log(f"pairs {len(keys):,} | test genes {len(test_gi):,} | phi q01/50/99 {np.quantile(phi,[.01,.5,.99]).round(4)}")
S_l = burden(phi); S_f = burden(np.ones_like(phi)); log("burden done")
Zl = {t: zrow(S_l[test_gi][:, TR[t]["cols"]]) for t in TRAITS}
Zf = {t: zrow(S_f[test_gi][:, TR[t]["cols"]]) for t in TRAITS}
del S_l, S_f
def metric(Z, r): return float(np.abs(Z @ r).mean() / (len(r) * r.std() + 1e-8))
obs = {t: metric(Zl[t], TR[t]["r"]) / metric(Zf[t], TR[t]["r"]) for t in TRAITS}
obs["mean"] = float(np.mean([obs[t] for t in TRAITS]))
log("observed ratio (v3 phi, fold0 holdout):", {k: round(v, 4) for k, v in obs.items()})
rng = np.random.default_rng(20260906); null = {t: [] for t in TRAITS}; null["mean"] = []
for i in range(NPERM):
    row = {}
    for t in TRAITS:
        rp = rng.permutation(TR[t]["r"])
        row[t] = metric(Zl[t], rp) / metric(Zf[t], rp)
        null[t].append(row[t])
    null["mean"].append(float(np.mean(list(row.values()))))
    if (i + 1) % 10 == 0: log(f"perm {i+1}/{NPERM} mean-ratio running median {np.median(null['mean']):.4f}")
res = {}
for k in list(TRAITS) + ["mean"]:
    a = np.array(null[k])
    res[k] = dict(obs=obs[k], null_mean=float(a.mean()), null_sd=float(a.std()),
                  q05=float(np.quantile(a, .05)), q50=float(np.quantile(a, .5)), q95=float(np.quantile(a, .95)),
                  mde_alpha05=float(np.quantile(a, .95)), z_obs=float((obs[k] - a.mean()) / (a.std() + 1e-12)),
                  p_perm=float((np.sum(a >= obs[k]) + 1) / (len(a) + 1)))
out = dict(nperm=NPERM, fold=0, phi="spline_fold0.v3", n_test_genes=int(len(test_gi)), pairs=int(len(keys)), result=res)
json.dump(out, open(_config_path("${PROJECT_ROOT}/work/ref/annot/l1/mde_fold0.json"), "w"), indent=1)
log("RESULT", json.dumps({k: {kk: round(vv, 4) for kk, vv in v.items()} for k, v in res.items()}))
print("MDE_DONE")
