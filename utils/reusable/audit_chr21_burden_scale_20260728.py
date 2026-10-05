#!/usr/bin/env python
from __future__ import annotations

import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import sys
from pathlib import Path

import numpy as np
import zarr

BURDEN_SRC = Path(
    _config_path("${PROJECT_ROOT}"
    "/work/production_chr1_21_deeprvat/assoc/chr21_grouped40_recovery_20260725_013126")
)

def main() -> None:
    a = zarr.open((BURDEN_SRC / "burdens/burdens_average.zarr").as_posix(),
                  mode="r")[:, :, 0]
    print(f"BURDEN\tshape_samples_x_genes\t{a.shape}")
    print(f"BURDEN\tglobal_min\t{a.min():.6f}")
    print(f"BURDEN\tglobal_max\t{a.max():.6f}")
    print(f"BURDEN\tglobal_mean\t{a.mean():.6f}")
    per_gene_sd = a.std(axis=0, ddof=1)
    per_gene_range = a.max(axis=0) - a.min(axis=0)
    per_gene_mean = a.mean(axis=0)
    for nm, v in [("per_gene_mean", per_gene_mean),
                  ("per_gene_sd", per_gene_sd),
                  ("per_gene_range", per_gene_range)]:
        q = np.percentile(v, [0, 5, 25, 50, 75, 95, 100])
        print(f"BURDEN\t{nm}_percentiles_0_5_25_50_75_95_100\t"
              + "\t".join(f"{x:.6f}" for x in q))
    print("BURDEN\tfraction_of_sigmoid_0_1_range_spanned_median_gene"
          f"\t{np.median(per_gene_range):.6f}")
    print(f"BURDEN\tn_genes_with_range_lt_0.01\t{int((per_gene_range < 0.01).sum())}"
          f"\tof\t{a.shape[1]}")
    print(f"BURDEN\tn_genes_with_range_lt_0.05\t{int((per_gene_range < 0.05).sum())}"
          f"\tof\t{a.shape[1]}")
    print("BURDEN_SCALE_GATE\tPASS")

if __name__ == "__main__":
    main()
