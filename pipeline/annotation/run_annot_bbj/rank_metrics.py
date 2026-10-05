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


import sys, os, json, time, hashlib
MODEL_DIR = _config_path("${PROJECT_ROOT}/work/run_band15/model_v10_out"); E6 = f"{MODEL_DIR}/e6_runs"
sys.path.insert(0, MODEL_DIR)
import numpy as np
from scipy.special import expit, logit
from scipy.optimize import minimize
from scipy.stats import rankdata
import l1_train_v10 as M
RUN = sys.argv[1]; extra = sys.argv[2:]
OUT = f"{E6}/{RUN}"; T0 = time.time()
def log(m): print(f"[rank {time.time()-T0:6.0f}s] {m}", flush=True)
L1 = _config_path("${PROJECT_ROOT}/work/ref/annot/bbj_l1")
cli = ["--out", f"{OUT}/rank_tmp", "--bbj-annot", f"{L1}/bbj_annot.tsv.gz", "--bbj-pip", f"{L1}/bbj_pip", "--traits", f"{L1}/traits.tsv",
       "--annot-stats", _config_path("${PROJECT_ROOT}/work/ref/annot/cache/fm_all.stats.json"), "--threads", "4", "--memory-gb", "60", "--tmpdir", f"{MODEL_DIR}/.scratch_e6/rank_{RUN}",
       "--allow-concurrent-heavy", "--c1-columns", "v9-eight", "--unseen-intercept-rule", "mean-trained", "--export-intercept-rule", "mean-trained",
       "--alpha", "0.05", "--power-target", "0.8", "--type1-max", "0.10"] + (extra or ["--fold", "0", "--phi-columns", "v9-eighteen", "--c1-ecdf-reference", "training", "--train-bands", "all"])
args = M.parser().parse_args(cli); args.custom_grid = False
os.makedirs(args.tmpdir, exist_ok=True)
M.load_libraries(args.threads)
res = json.load(open(f"{OUT}/RESULT.json")); model = json.load(open(f"{OUT}/phi.json")); fp = np.load(f"{OUT}/fitted_parameters.npz")
data = M.Data(args); M._trim_memory(); log(f"data loaded: variants {len(data.keys):,} rows {len(data.y):,}")
splits = M.split_data(args, data); train, idx = splits["train"], splits["chromosome"]
for k in ("train", "chromosome"):
    assert hashlib.sha256(splits[k].tobytes()).hexdigest() == res["split_hashes"][k], f"GATE FAIL: split hash {k}"
log(f"splits verified: train {len(train):,} test {len(idx):,}")
y = data.y; trait_ids = sorted(set(data.ti[train])); tname = {i: t["trait"] for i, t in enumerate(data.traits)}
basis = M.transform_design(model["design"], data.X, data.cols); coef = np.asarray(model["coefficients"]["shared"])
assert basis.shape[1] == len(coef), "GATE FAIL: basis/coef width"
def phi_pred(indices):
    eta = np.empty(len(indices)); step = 500000
    for lo in range(0, len(indices), step):
        ii = indices[lo:lo+step]; eta[lo:lo+step] = basis[data.vi[ii]].astype(float) @ coef
    ti = data.ti[indices]; inter = np.full(len(indices), float(model.get("global_intercept", np.mean(list(model["trait_intercepts"].values())))))
    for i, name in tname.items():
        if name in model["trait_intercepts"]: inter[ti == i] = model["trait_intercepts"][name]
    return M.open_probability(expit(eta + inter))
c2_design = M.make_design(args, data, train, True); c2_a = fp["C2"]; ncoef2 = c2_design.nparam
assert len(c2_a) == ncoef2 + len(trait_ids), f"GATE FAIL: C2 param length {len(c2_a)} vs {ncoef2}+{len(trait_ids)}"
def c2_pred(indices):
    eta = c2_design.B[indices].astype(float) @ c2_a[:ncoef2]; ti = data.ti[indices]; unknown = np.ones(len(indices), bool)
    for j, t in enumerate(trait_ids):
        m = ti == t; eta[m] += c2_a[ncoef2 + j]; unknown[m] = False
    if unknown.any(): eta[unknown] += float(np.mean(c2_a[ncoef2:]))
    return M.open_probability(expit(eta))
P = {}
P["phi"] = phi_pred(idx); P["C2"] = c2_pred(idx)
c1, _ = M.c1_predictions(args, data, train, idx); P["C1"] = c1["C1"]
cal = M.c1_calibrate(args, data, train, y); P["C1_cal"] = M.apply_c1_calibration(cal, P["C1"])
rc = res["evaluation"]["chromosome"]["C1_calibration"]
log(f"C1_cal refit a={cal['intercept']:.6f} b={cal['slope']:.6f} vs RESULT a={rc['intercept']:.6f} b={rc['slope']:.6f}")
assert abs(cal["intercept"] - rc["intercept"]) < 1e-6 and abs(cal["slope"] - rc["slope"]) < 1e-6, "GATE FAIL: C1 calibration mismatch (prediction pipeline not reproduced)"
def wce(p, yy, w):
    return float(np.sum(w * (-(yy * np.log(p) + (1 - yy) * np.log(1 - p)))))
bands_idx = data.bands[idx]; checks = {}
for b in M.BANDS:
    m = bands_idx == b; ib = idx[m]; w, _ = M.weights(data, ib, y)
    checks[b] = {k: (wce(P[k][m], y[ib], w), res["evaluation"]["chromosome"]["bands"][b]["models"][k]["weighted_ce"]["observed"]) for k in ("phi", "C2", "C1", "C1_cal")}
log("CE check (recomputed, RESULT): " + json.dumps({b: {k: [round(v[0], 5), round(v[1], 5)] for k, v in d.items()} for b, d in checks.items()}))
ce_ok = all(abs(v[0] - v[1]) < 1e-3 for d in checks.values() for v in d.values())
log("stack fit ...")
ptr_c2 = c2_pred(train); ptr_phi = phi_pred(train)
Xs = np.column_stack([np.ones(len(train)), logit(np.clip(ptr_c2, M.NUM_EPS, 1 - M.NUM_EPS)), logit(np.clip(ptr_phi, M.NUM_EPS, 1 - M.NUM_EPS))]); del ptr_c2, ptr_phi
wtr, _ = M.weights(data, train, y); ytr = y[train]
def obj(th):
    eta = Xs @ th; loss = np.sum(wtr * (np.logaddexp(0., eta) - ytr * eta)); r = wtr * (expit(eta) - ytr); return float(loss), Xs.T @ r
st = minimize(obj, np.array([logit(np.clip(np.average(ytr, weights=wtr), 1e-6, 1-1e-6)), 1., 0.]), jac=True, method="L-BFGS-B", options=dict(maxiter=500, ftol=1e-12, gtol=1e-10))
stack = dict(coef=[float(v) for v in st.x], success=bool(st.success), nit=int(st.nit), features=["1", "logit(C2)", "logit(phi)"], n_train=int(len(train)))
log(f"stack {stack}"); del Xs, wtr
P["stack"] = M.open_probability(expit(st.x[0] + st.x[1] * logit(np.clip(P["C2"], M.NUM_EPS, 1 - M.NUM_EPS)) + st.x[2] * logit(np.clip(P["phi"], M.NUM_EPS, 1 - M.NUM_EPS))))
MODELS = ["phi", "C1", "C1_cal", "C2", "stack"]
B_PERM = 50; rng = np.random.default_rng(123)
def spearman_groups(pred, lab, groups):
    out = {}
    for g, ii in groups.items():
        if len(ii) < 3: continue
        l = lab[ii]; p = pred[ii]
        if np.ptp(l) == 0 or np.ptp(p) == 0: continue
        rl = rankdata(l); rp = rankdata(p); rl -= rl.mean(); rp -= rp.mean()
        out[g] = float(np.dot(rl, rp) / np.sqrt(np.dot(rl, rl) * np.dot(rp, rp)))
    return out
def agg(rhos, groups_trait):
    if not rhos: return dict(trait_mean=None, region_mean=None, n_regions=0, n_traits=0)
    per_trait = {}
    for g, r in rhos.items(): per_trait.setdefault(groups_trait[g], []).append(r)
    tm = float(np.mean([np.mean(v) for v in per_trait.values()]))
    return dict(trait_mean=tm, region_mean=float(np.mean(list(rhos.values()))), n_regions=len(rhos), n_traits=len(per_trait))
rows = []; detail = {}
ti_idx = data.ti[idx]; reg_idx = data.regions[idx]
for b in list(M.BANDS) + ["all"]:
    m = np.ones(len(idx), bool) if b == "all" else (bands_idx == b)
    sel = np.where(m)[0]
    groups = {}
    for k in sel: groups.setdefault((int(ti_idx[k]), reg_idx[k]), []).append(k)
    groups = {g: np.asarray(v) for g, v in groups.items()}; gtrait = {g: g[0] for g in groups}
    lab = y[idx]
    obs = {mn: agg(spearman_groups(P[mn], lab, groups), gtrait) for mn in MODELS}
    null = {mn: [] for mn in MODELS}; null_d = []
    for r in range(B_PERM):
        lp = lab.copy()
        for g, ii in groups.items():
            if len(ii) > 1: lp[ii] = lab[ii][rng.permutation(len(ii))]
        a = {mn: agg(spearman_groups(P[mn], lp, groups), gtrait)["trait_mean"] for mn in MODELS}
        for mn in MODELS: null[mn].append(a[mn] if a[mn] is not None else np.nan)
        null_d.append((a["stack"] - a["C2"]) if (a["stack"] is not None and a["C2"] is not None) else np.nan)
    for mn in MODELS:
        o = obs[mn]; nm = np.nanmean(null[mn]) if o["n_regions"] else None; ns = np.nanstd(null[mn], ddof=1) if o["n_regions"] else None
        rows.append(dict(run=RUN, band=b, metric="region_spearman_trait_mean", model=mn, value=o["trait_mean"], n=o["n_regions"], n_traits=o["n_traits"],
                         null_mean=nm, null_sd=ns, z=((o["trait_mean"] - nm) / ns if (o["trait_mean"] is not None and ns) else None), extra=dict(region_mean=o["region_mean"])))
    dobs = (obs["stack"]["trait_mean"] - obs["C2"]["trait_mean"]) if (obs["stack"]["trait_mean"] is not None and obs["C2"]["trait_mean"] is not None) else None
    nm = np.nanmean(null_d) if dobs is not None else None; ns = np.nanstd(null_d, ddof=1) if dobs is not None else None
    rows.append(dict(run=RUN, band=b, metric="delta_spearman_stack_minus_C2", model="stack-C2", value=dobs, n=obs["stack"]["n_regions"], n_traits=obs["stack"]["n_traits"], null_mean=nm, null_sd=ns, z=((dobs - nm) / ns if (dobs is not None and ns) else None), extra={}))
    pos = lab[sel] >= 0.1; npos = int(pos.sum())
    for mn in MODELS:
        if npos == 0 or npos == len(sel): ap = None
        else:
            order = np.argsort(-P[mn][sel], kind="mergesort"); tp = np.cumsum(pos[order]); prec = tp / np.arange(1, len(sel) + 1); rec = tp / npos
            ap = float(np.sum(np.diff(np.concatenate([[0.], rec])) * prec))
        rows.append(dict(run=RUN, band=b, metric="auprc_pip_ge_0.1", model=mn, value=ap, n=int(len(sel)), n_traits=npos, null_mean=(npos / len(sel) if len(sel) else None), null_sd=None, z=None, extra=dict(baseline="positive_rate")))
    if b != "all":
        for mn in ("phi", "C1", "C1_cal", "C2"):
            s = res["evaluation"]["chromosome"]["bands"][b]["models"][mn]["cs_spearman"]
            rows.append(dict(run=RUN, band=b, metric="cs_spearman_old", model=mn, value=s.get("observed"), n=res["evaluation"]["chromosome"]["bands"][b]["models"][mn].get("cs_defined"), n_traits=None, null_mean=s.get("null_mean"), null_sd=s.get("null_sd"), z=s.get("z"), extra={}))
    log(f"band {b}: rows {len(sel):,} regions {len(groups)} spearman(trait_mean) " + json.dumps({mn: (None if obs[mn]['trait_mean'] is None else round(obs[mn]['trait_mean'], 4)) for mn in MODELS}))
keys = ["run", "band", "metric", "model", "value", "n", "n_traits", "null_mean", "null_sd", "z", "extra"]
with open(f"{OUT}/summary_rank.tsv.tmp", "w") as o:
    o.write("\t".join(keys) + "\n")
    for r in rows: o.write("\t".join(json.dumps(r[k]) if isinstance(r[k], dict) else ("" if r[k] is None else (f"{r[k]:.6g}" if isinstance(r[k], float) else str(r[k]))) for k in keys) + "\n")
os.replace(f"{OUT}/summary_rank.tsv.tmp", f"{OUT}/summary_rank.tsv")
json.dump(dict(run=RUN, ce_check=checks, ce_check_ok=ce_ok, c1_calibration=cal, stack=stack, permutations=B_PERM, seed=123,
               null_definition="labels permuted within (trait, region) blocks among band rows", spearman_rule="(trait,region) groups n>=3, non-constant; trait_mean = mean over traits of mean over that trait's regions",
               prereg_note="D1' rank metrics defined after inspecting primary results (v83); stacking and null are this script's definitions (STAR)"),
          open(f"{OUT}/summary_rank.meta.json", "w"), indent=1)
assert ce_ok, "GATE FAIL: recomputed CE differs from RESULT by >1e-3 (see meta)"
open(f"{OUT}/RANK_DONE", "w").write("ok\n"); log("RANK_DONE")
