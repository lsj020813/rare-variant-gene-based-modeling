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


import shutil, sys, hashlib, py_compile
P = _config_path("${PROJECT_ROOT}/work/run_band15/model_v10_out/l1_train_v10.py")
src = open(P).read()
if "E6D_PATCH" in src: print("already patched"); sys.exit(0)
assert "E6C_PATCH" in src
shutil.copy(P, P + ".pre_e6d_bak")
def rep(old, new):
    global src
    assert src.count(old) == 1, f"count {src.count(old)}: {old[:70]}"
    src = src.replace(old, new)
rep("class PermutationPlan:\n",
    "def c1_calibrate(args, data, train, y):\n"
    "    \"\"\"E6D_PATCH: two-parameter logistic recalibration of the C1 Cauchy score on training rows (same weights as the arms).\"\"\"\n"
    "    c1, _ = c1_predictions(args, data, train, train)\n"
    "    x = logit(np.clip(c1['C1'], NUM_EPS, 1 - NUM_EPS))\n"
    "    w, _ = weights(data, train, y)\n"
    "    yt = y[train]\n"
    "    def obj(theta):\n"
    "        eta = theta[0] + theta[1] * x\n"
    "        loss = np.sum(w * (np.logaddexp(0., eta) - yt * eta))\n"
    "        r = w * (expit(eta) - yt)\n"
    "        return float(loss), np.array([r.sum(), (r * x).sum()])\n"
    "    start = np.array([logit(np.clip(np.average(yt, weights=w), NUM_EPS, 1 - NUM_EPS)), 0.])\n"
    "    res = minimize(obj, start, jac=True, method='L-BFGS-B', options=dict(maxiter=500, ftol=1e-12, gtol=1e-10))\n"
    "    require(np.isfinite(res.x).all(), 'C1 calibration failed')\n"
    "    return dict(intercept=float(res.x[0]), slope=float(res.x[1]), success=bool(res.success), nit=int(res.nit), n_train=int(len(train)),\n"
    "                fitted_on='outer training rows (same rows as final C2/phi fits)', model='expit(a + b*logit(C1_cauchy))')\n\n\n"
    "def apply_c1_calibration(cal, c1_scores):\n"
    "    return open_probability(expit(cal['intercept'] + cal['slope'] * logit(np.clip(c1_scores, NUM_EPS, 1 - NUM_EPS))))\n\n\n"
    "class PermutationPlan:\n")
rep("    flat = float(np.mean(y[train]))\n    reports = {}\n",
    "    flat = float(np.mean(y[train]))\n    c1_cal = c1_calibrate(args, data, train, y)  # E6D_PATCH\n    reports = {}\n")
rep("        predictions = dict(flat=np.full(len(idx), flat), C2=final_c2.predict(c2_a, idx),\n                           **predictions, phi=final_phi.predict(phi_a, idx))\n        reports[axis], _ = evaluate(args, data, idx, predictions, args.test_seed, y)\n        reports[axis]['C1_specification'] = c1_meta\n",
    "        predictions = dict(flat=np.full(len(idx), flat), C2=final_c2.predict(c2_a, idx),\n                           **predictions, C1_cal=apply_c1_calibration(c1_cal, predictions['C1']), phi=final_phi.predict(phi_a, idx))  # E6D_PATCH\n        reports[axis], _ = evaluate(args, data, idx, predictions, args.test_seed, y)\n        reports[axis]['C1_specification'] = c1_meta\n        reports[axis]['C1_calibration'] = c1_cal  # E6D_PATCH\n")
rep("        for left, right in zip(sequence, sequence[1:]):\n            diff = ces[:, names.index(left)] - ces[:, names.index(right)]\n            ladder[left + '->' + right] = null_summary(float(diff[0]), diff[1:])\n",
    "        if 'C1_cal' in names and 'C1' in names and 'phi' in names:  # E6D_PATCH: recalibrated C1 steps alongside the original ladder\n"
    "            sequence = sequence + [('C1', 'C1_cal'), ('C1_cal', 'phi')]\n"
    "        pairs = [p for p in zip(sequence, sequence[1:]) if isinstance(p[0], str) and isinstance(p[1], str)] + [p for p in sequence if isinstance(p, tuple)]\n"
    "        for left, right in pairs:\n            diff = ces[:, names.index(left)] - ces[:, names.index(right)]\n            ladder[left + '->' + right] = null_summary(float(diff[0]), diff[1:])\n")
open(P, "w").write(src); py_compile.compile(P, doraise=True)
print("PATCHED4", hashlib.sha256(open(P,'rb').read()).hexdigest()[:16], src.count("E6D_PATCH"))
