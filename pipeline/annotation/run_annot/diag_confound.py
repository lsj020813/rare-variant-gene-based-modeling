
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import numpy as np, json, sys, os, re
sys.argv = ["x", "--smoke", "22", "--arm", "spline", "--gradcheck-only"]
src = open("l1_train.py").read()
cut = src.index("if A.gcv_probe:")
g = {}; exec(compile(src[:cut], "l1_train_head", "exec"), g)
B, X, cols, TR, TRAITS = g["B"], g["X"], g["cols"], g["TR"], g["TRAITS"]
test_gi, train_gi, GID, burden, zrow = g["test_gi"], g["train_gi"], g["GID"], g["burden"], g["zrow"]
a = np.load(_config_path("${PROJECT_ROOT}/work/ref/annot/l1/spline_smoke_chr22.a.npy")); phi = np.load(_config_path("${PROJECT_ROOT}/work/ref/annot/l1/spline_smoke_chr22.phi.npy"))
S_l = burden(phi.astype(np.float32)); S_f = burden(np.ones(len(phi), np.float32))
out = {}
for t in TRAITS:
    c, r = TR[t]["cols"], TR[t]["r"]; n = len(c); sig = r.std()
    Tl = (zrow(S_l[test_gi][:, c]) @ r) / (np.sqrt(n) * sig); Tf = (zrow(S_f[test_gi][:, c]) @ r) / (np.sqrt(n) * sig)
    d = np.abs(Tl) - np.abs(Tf); o = np.argsort(-d)
    out[t] = dict(ratio=float(np.abs(Tl).mean() / np.abs(Tf).mean()), medT2_flat=float(np.median(Tf**2)), medT2_learn=float(np.median(Tl**2)),
                  top5_share=float(d[o[:5]].sum() / max(d.sum(), 1e-9)), n_up=int((d > 0).sum()), n_test=len(d),
                  top5=[(int(test_gi[i]), round(float(Tf[i]), 2), round(float(Tl[i]), 2)) for i in o[:5]])
cont = [j for j, cn in enumerate(cols) if cn not in g["BINARY"]]; binc = [j for j, cn in enumerate(cols) if cn in g["BINARY"]]
load = {}
k = 0
for j in cont: load[cols[j]] = float(np.abs(a[k:k+7]).sum()); k += 7
for j in binc: load[cols[j]] = float(abs(a[k])); k += 1
top = sorted(load.items(), key=lambda kv: -kv[1])[:12]
tech = [cn for cn in cols if cn in ("maf", "r2", "avg_cs", "is_typed", "is_indel")]
corr = {cn: float(np.corrcoef(np.nan_to_num(X[:, cols.index(cn)]), np.log(phi))[0, 1]) for cn in tech}
print(json.dumps(dict(per_trait=out, top_loadings=top, corr_logphi_tech=corr, phi_q=[float(q) for q in np.quantile(phi, [0, .01, .5, .99, 1])]), indent=1))
