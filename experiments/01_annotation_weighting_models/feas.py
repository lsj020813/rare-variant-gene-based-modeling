import os as _cfg_os
import math as _cfg_math

def _config_number(name, cast=int, positive=False):
    value = cast(_cfg_os.environ[name])
    if not _cfg_math.isfinite(value) or (positive and value <= 0):
        raise ValueError(name + " must be finite" + (" and positive" if positive else ""))
    return value

import os as _os
N_SAMPLES = _config_number("N_SAMPLES", int, True)

import torch
V, N = 36957, N_SAMPLES
for dt, name in ((torch.float32,"fp32"), (torch.float16,"fp16")):
    bytes_D = V * N * torch.tensor([], dtype=dt).element_size()
    print(f"{name}: D = {bytes_D/1e9:.1f} GB  (A2 capacity 14.6 GB)")
print("stream design: D stays on CPU RAM (float16 = %.1f GB, server RAM 250GB ok)" % (V*N*2/1e9))
print("GPU holds only per-block slices + S matrix (6959 x N_SAMPLES fp32 = %.1f GB)" % (6959*N_SAMPLES*4/1e9))
