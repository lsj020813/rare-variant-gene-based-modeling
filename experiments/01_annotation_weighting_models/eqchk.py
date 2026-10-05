
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import numpy as np
v1 = np.load(_config_path("${PROJECT_ROOT}/work/ref/l1/B/chr21.B.npz"))
v2 = np.load(_config_path("${PROJECT_ROOT}/work/run_l1/v2check/chr21.B.npz"))
r1, r2 = v1["rows"], v2["rows"]
print("rows identical:", r1.shape == r2.shape and bool((r1 == r2).all()), "| n =", r1.shape[0])
B1, B2 = v1["B"].astype(np.float32), v2["B"].astype(np.float32)
print("shape:", B1.shape, B2.shape)
d = np.abs(B1 - B2)
print("max abs diff:", float(d.max()), "| mean:", float(d.mean()), "| exact equal cells:", float((d == 0).mean()))
print("EQUIV_OK" if d.max() <= 1e-3 else "EQUIV_FAIL")
