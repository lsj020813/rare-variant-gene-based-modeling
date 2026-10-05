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
if "E6F_PATCH" in src: print("already patched"); sys.exit(0)
assert "E6E_PATCH" in src
shutil.copy(P, P + ".pre_e6f_bak")
def rep(old, new):
    global src
    assert src.count(old) == 1, f"count {src.count(old)}: {old[:70]}"
    src = src.replace(old, new)
rep("def weights(data, indices, y=None):\n", "def weights(data, indices, y=None, base=None):\n")
rep("    base = np.ones(len(indices), dtype=float)\n    blocks = make_blocks([(int(data.ti[i]), data.regions[i], int(data.cs[i])) for i in indices], data.cs[indices] != -1)\n    for block in blocks:\n        base[block] = 1 / len(block)\n",
    "    if base is None:  # E6F_PATCH: label-independent block mass may be supplied precomputed\n"
    "        base = weights_base(data, indices)\n    base = np.asarray(base, dtype=float)\n")
rep("def split_data(args, data):\n",
    "def weights_base(data, indices):\n"
    "    \"\"\"E6F_PATCH: CS mass=1 / non-CS mass=1 base weights (independent of labels).\"\"\"\n"
    "    indices = np.asarray(indices, dtype=int)\n"
    "    base = np.ones(len(indices), dtype=float)\n"
    "    blocks = make_blocks([(int(data.ti[i]), data.regions[i], int(data.cs[i])) for i in indices], data.cs[indices] != -1)\n"
    "    for block in blocks:\n        base[block] = 1 / len(block)\n"
    "    return base\n\n\ndef split_data(args, data):\n")
rep("        ces, ranks, counts = [], [], []\n        for yy in ys:\n            # Derive weights from this label vector without constructing outside-region labels.\n            work_y = y.copy()\n            work_y[indices] = yy\n            w, wm = weights(data, idx, work_y)\n",
    "        ces, ranks, counts = [], [], []\n        base_w = weights_base(data, idx)  # E6F_PATCH: computed once per band\n"
    "        work_y = y.copy()\n"
    "        for yy in ys:\n            # Derive weights from this label vector without constructing outside-region labels.\n            work_y[indices] = yy\n            w, wm = weights(data, idx, work_y, base_w)\n")
rep("        observed_w, observed_wm = weights(data, idx, y)\n", "        observed_w, observed_wm = weights(data, idx, y, base_w)\n")
rep("    ys = [y[indices]] + list(plan.labels(y))\n",
    "    def label_iter():  # E6F_PATCH: regenerate the deterministic permutation stream per band instead of holding B+1 label vectors\n"
    "        yield y[indices]\n        yield from plan.labels(y)\n")
rep("        for yy in ys:\n            # Derive weights", "        for yy in label_iter():\n            # Derive weights")
open(P, "w").write(src); py_compile.compile(P, doraise=True)
print("PATCHED6", hashlib.sha256(open(P,'rb').read()).hexdigest()[:16], src.count("E6F_PATCH"))
