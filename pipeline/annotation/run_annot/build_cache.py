#!/usr/bin/env python3
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

import os, sys, json, numpy as np
R = _config_path("${PROJECT_ROOT}/work/ref/annot"); OUT = f"{R}/cache"; os.makedirs(OUT, exist_ok=True)
FOLD = {7:0,13:0,16:0,19:0,21:0, 2:1,8:1,9:1,10:1,15:1, 1:2,6:2,14:2, 3:3,12:3,17:3,22:3, 4:4,5:4,11:4,18:4,20:4}
def require(c, m):
    if not c: raise RuntimeError(f"GATE FAIL: {m}")
keys, genes, folds, chrs, blocks = [], [], [], [], []
NONUNI = []
cols = None
for N in range(1, 23):
    fm = f"{R}/fm/chr{N}.fm.tsv"; t4 = f"{R}/t4/chr{N}.t4.tsv"
    require(os.path.exists(fm + ".done") and os.path.exists(t4 + ".done"), f"chr{N} inputs incomplete")
    T4 = {}
    with open(t4) as f:
        h4 = next(f).rstrip("\n").split("\t"); require(h4 == ["key37","maf","r2","avg_cs","is_typed","is_indel"], f"t4 header chr{N}")
        for line in f:
            a = line.rstrip("\n").split("\t"); T4[a[0]] = a[1:]
    with open(fm) as f:
        hdr = next(f).rstrip("\n").split("\t")
        na_cols = [c for c in hdr if c.endswith("_na")]
        val_cols = [c for c in hdr if c not in ("key37","gene","gene_base") and not c.endswith("_na")]
        newcols = val_cols + ["t1_na"] + h4[1:]
        if cols is None: cols = newcols
        require(cols == newcols, f"schema drift chr{N}")
        ix = {c: i for i, c in enumerate(hdr)}
        rows = []
        miss_t4 = 0; nonuni = 0
        for line in f:
            a = line.rstrip("\n").split("\t")
            na = [a[ix[c]] for c in na_cols]
            if len(set(na)) != 1: nonuni += 1
            na = ["1" if "1" in na else "0"]
            t4v = T4.get(a[0])
            if t4v is None: miss_t4 += 1; t4v = [""] * 5
            vals = [a[ix[c]] for c in val_cols] + [na[0]] + t4v
            rows.append([float(x) if x != "" else np.nan for x in vals])
            keys.append(a[0]); genes.append(a[ix["gene_base"]]); folds.append(FOLD[N]); chrs.append(N)
        require(miss_t4 == 0, f"chr{N}: {miss_t4} rows lack T4 (band VCF should cover all assigned variants)")
        require(nonuni <= 100, f"chr{N}: {nonuni} rows with non-uniform T1 mask (expected ~16 genome-wide)")
        NONUNI.append(nonuni)
        blocks.append(np.asarray(rows, dtype=np.float32))
    print(f"[chr{N}] {len(rows)} rows", flush=True)
X = np.vstack(blocks)
require(X.shape[0] == _config_number("N_ANNOTATION_VARIANTS", int, True), f"row total {X.shape[0]} != configured expected row count")
np.savez(f"{OUT}/fm_all.npz", X=X, cols=np.array(cols), key37=np.array(keys), gene=np.array(genes),
         fold=np.array(folds, dtype=np.int8), chr=np.array(chrs, dtype=np.int8))
na_pct = {c: round(float(np.isnan(X[:, i]).mean() * 100), 2) for i, c in enumerate(cols)}
json.dump({"shape": X.shape, "cols": cols, "na_pct": na_pct, "t1_mask_nonuniform_rows": int(sum(NONUNI))}, open(f"{OUT}/fm_all.stats.json", "w"))
print("CACHE_DONE", X.shape, f"{X.nbytes/1e9:.2f} GB")
