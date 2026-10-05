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
if "E6G_PATCH" in src: print("already patched"); sys.exit(0)
assert "E6F_PATCH" in src
shutil.copy(P, P + ".pre_e6g_bak")
def rep(old, new):
    global src
    assert src.count(old) == 1, f"count {src.count(old)}: {old[:70]}"
    src = src.replace(old, new)
rep("    ap.add_argument('--inner-va-folds', type=int, choices=[1, 2], default=2)\n",
    "    ap.add_argument('--inner-va-folds', type=int, choices=[1, 2], default=2)\n"
    "    ap.add_argument('--folds-override', help='E6G_PATCH: JSON {fold: [chromosomes]} replacing FOLDS (subset simulations); SEALED unchanged')\n")
rep("    args = parser().parse_args()\n    args.custom_grid = '--lam-divisors' in sys.argv\n",
    "    args = parser().parse_args()\n    args.custom_grid = '--lam-divisors' in sys.argv\n"
    "    if getattr(args, 'folds_override', None):  # E6G_PATCH\n"
    "        override = json.loads(args.folds_override)\n"
    "        require(sorted(int(k) for k in override) == [0, 1, 2, 3, 4], 'folds-override must define folds 0..4')\n"
    "        FOLDS.clear(); FOLDS.update({int(k): set(str(c) for c in v) for k, v in override.items()})\n"
    "        print(f'[E6G_PATCH] FOLDS overridden: {sorted((k, sorted(v)) for k, v in FOLDS.items())}', flush=True)\n")
rep("    result = dict(version='v10', run_kind='external_PIP_annotation_only', fold=args.fold,\n",
    "    result = dict(version='v10', run_kind='external_PIP_annotation_only', fold=args.fold, folds_used={str(k): sorted(v) for k, v in FOLDS.items()},  # E6G_PATCH\n")
open(P, "w").write(src); py_compile.compile(P, doraise=True)
print("PATCHED7", hashlib.sha256(open(P,'rb').read()).hexdigest()[:16], src.count("E6G_PATCH"))
