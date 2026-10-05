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


import re, shutil, sys, hashlib
P = _config_path("${PROJECT_ROOT}/work/run_band15/model_v10_out/l1_train_v10.py")
src = open(P).read()
if "E6_PATCH" in src: print("already patched"); sys.exit(0)
shutil.copy(P, P + ".pre_e6bak")
def rep(old, new):
    global src
    assert src.count(old) == 1, f"pattern count {src.count(old)}: {old[:80]}"
    src = src.replace(old, new)
rep("    ap.add_argument('--phi-columns', choices=['all34', 'v9-eighteen'], help='★ 34-column brief vs v9 18-column design; no default')",
    "    ap.add_argument('--phi-columns', help=\"★ 'all34' | 'v9-eighteen' (S1+cadd+S2+S3=18) | 'v8-seventeen' (S1+S2+S3, no cadd) | 'list:c1,c2,...' explicit; no default  [E6_PATCH 2026-09-08]\")")
rep("        selected = self.cols if control or feature_mode == 'all34' else S1 + ['cadd'] + S2 + S3\n",
    "        selected = phi_feature_set(self.cols, feature_mode, control)  # E6_PATCH\n")
rep("class Design:\n",
    "def phi_feature_set(cols, feature_mode, control):\n"
    "    \"\"\"E6_PATCH 2026-09-08: explicit phi column sets. v9-eighteen already contains cadd (v9 appended cadd_side to the 17-column v8 set).\"\"\"\n"
    "    if control or feature_mode == 'all34':\n        return list(cols)\n"
    "    if feature_mode == 'v9-eighteen':\n        return S1 + ['cadd'] + S2 + S3\n"
    "    if feature_mode == 'v8-seventeen':\n        return S1 + S2 + S3\n"
    "    require(isinstance(feature_mode, str) and feature_mode.startswith('list:'), 'unknown --phi-columns mode')\n"
    "    chosen = [c for c in feature_mode[5:].split(',') if c]\n"
    "    require(chosen and len(set(chosen)) == len(chosen) and set(chosen) <= set(cols), 'explicit phi column list must be unique and within the 34 annotation columns')\n"
    "    return chosen\n\n\nclass Design:\n")
rep("        n = integer(row['N'], 'trait N')\n        require(n > 0, 'trait N must be positive')\n",
    "        if str(row['N']).strip().upper() in ('NA', ''):  # E6_PATCH: N unavailable from BBJ files -> metadata None\n            n = None\n        else:\n            n = integer(row['N'], 'trait N')\n            require(n > 0, 'trait N must be positive')\n")
rep("def oracle_signal(data, spec):\n    require(set(spec) >= {'terms', 'intercepts'} and spec['terms'], 'oracle specification schema')\n    signal = np.zeros(len(data.keys))\n",
    "def oracle_signal(data, spec):\n"
    "    if spec.get('phi_model'):  # E6_PATCH 2026-09-08: oracle phi* = exported fitted phi (design basis @ shared coefficients), scaled by delta downstream\n"
    "        require('intercepts' in spec, 'oracle trait intercepts missing')\n"
    "        with open(input_path(spec['phi_model'])) as fh:\n            m = json.load(fh)\n"
    "        require(m.get('version') == 'v10' and list(m['coefficients']) == ['shared'], 'phi_model oracle needs a shared-phi v10 export')\n"
    "        basis = transform_design(m['design'], data.X, data.cols)\n"
    "        signal = basis @ np.asarray(m['coefficients']['shared'])\n"
    "        require(set(spec['intercepts']) >= {r['trait'] for r in data.traits}, 'oracle trait intercepts missing')\n"
    "        return signal\n"
    "    require(set(spec) >= {'terms', 'intercepts'} and spec['terms'], 'oracle specification schema')\n    signal = np.zeros(len(data.keys))\n")
open(P, "w").write(src)
import py_compile; py_compile.compile(P, doraise=True)
print("PATCHED", hashlib.sha256(open(P,'rb').read()).hexdigest()[:16])
