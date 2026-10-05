import os as _cfg_os
import re as _cfg_re
import math as _cfg_math

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)

def _config_number(name, cast=int, positive=False):
    value = cast(_cfg_os.environ[name])
    if not _cfg_math.isfinite(value) or (positive and value <= 0):
        raise ValueError(name + " must be finite" + (" and positive" if positive else ""))
    return value

import os as _os
N_SAMPLES = _config_number("N_SAMPLES", int, True)

import numpy as np, sys, json, glob, os
sys.argv = ["l1_train_v7.py", "--fold", "0", "--arm", "spline"]
src = open(_config_path("${PROJECT_ROOT}/work/run_band15/l1_train_v7.py")).read()
import ast
nodes = ast.parse(src).body
start = next(n.lineno - 1 for n in nodes if isinstance(n, ast.Assign)
             and any(isinstance(t, ast.Name) and t.id == "DSC" for t in n.targets))
end = next(n.lineno - 1 for n in nodes if isinstance(n, ast.Assign)
           and any(isinstance(t, ast.Name) and t.id == "_ord" for t in n.targets))
lines = src.splitlines(keepends=True)
src = "".join(lines[:start]) + "NS = N_SAMPLES\n" + "".join(lines[end:])
cut = src.index("M_full = ")
g = {"__file__": _config_path("${PROJECT_ROOT}/work/run_band15/l1_train_v7.py")}
exec(compile(src[:cut], "trainer_prefix", "exec"), g)
B, X, cols, BLOCKS, CL, train_pairs, test_gi, GID, gl, gchr = (g[k] for k in ["B","X","cols","BLOCKS","CL","train_pairs","test_gi","GID","gl","gchr"])
cks = sorted(glob.glob(_config_path("${PROJECT_ROOT}/work/ref15/annot/l1/spline_v7_fold0*.ckpt.npz")), key=os.path.getmtime)
print("checkpoints:", [os.path.basename(c) for c in cks])
z = np.load(cks[-1], allow_pickle=True); a = z["a"].astype(np.float64); print("using", os.path.basename(cks[-1]), "nit", int(z["nit"]), "|a|max", round(float(np.abs(a).max()),3))
f = B.astype(np.float64) @ a; phi = np.exp(np.clip(f, -CL, CL))
out = dict(nit=int(z["nit"]), rms_f=float(np.sqrt((f[train_pairs]**2).mean())), clamp_frac_all=float((np.abs(f) >= CL).mean()), clamp_hi=float((f >= CL).mean()), clamp_lo=float((f <= -CL).mean()),
           phi_q=[float(q) for q in np.quantile(phi, [0,.01,.1,.5,.9,.99,1])], f_q=[float(q) for q in np.quantile(f, [0,.01,.5,.99,1])])
def colv(cn): return X[:, cols.index(cn)].astype(np.float64)
from scipy.stats import spearmanr
tech = {}
for cn in ["maf","r2","avg_cs","is_typed","is_indel","n_genes","gh_n_genes","dist_tss","cons","gpn_msa","linsight","epi_active","tf"]:
    if cn in cols:
        v = colv(cn); m = ~np.isnan(v); tech[cn] = dict(pearson_logphi=float(np.corrcoef(v[m], f[m])[0,1]), spearman=float(spearmanr(v[m], f[m]).correlation))
out["corr_f_vs_col"] = tech
maf = colv("maf"); qs = np.quantile(maf, [0,.25,.5,.75,1]); out["phi_mean_by_maf_quartile"] = [float(phi[(maf>=qs[i])&(maf<=qs[i+1])].mean()) for i in range(4)]
out["maf_quartile_edges"] = [float(q) for q in qs]
load = {}; k = 0
for b in BLOCKS:
    d = b[2].shape[1]; load[b[0]] = float(np.abs(a[k:k+d]).sum()); k += d
out["top_loadings"] = sorted(load.items(), key=lambda kv: -kv[1])[:12]
te = np.isin(GID, test_gi); out["clamp_frac_test_vs_train"] = [float((np.abs(f[te])>=CL).mean()), float((np.abs(f[~te])>=CL).mean())]
print("DIAG_JSON " + json.dumps(out))
