#!/usr/bin/env python3

from __future__ import annotations

import argparse
import gzip
import os
import sys

import numpy as np
import pandas as pd

MAF_MAX = 1e-3
CADD_MIN = 5.0
KEYS = ["chrom", "pos", "ref", "alt"]

def vcf_r2(path: str) -> pd.DataFrame:
    rows = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            info = f[7]
            r2 = None
            for kv in info.split(";"):
                if kv.startswith("R2="):
                    try:
                        r2 = float(kv[3:])
                    except ValueError:
                        r2 = None
                    break
            c = f[0][3:] if f[0].startswith("chr") else f[0]
            rows.append((int(c), int(f[1]), f[3], f[4], r2))
    return pd.DataFrame(rows, columns=KEYS + ["R2"])

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--annotations", required=True)
    ap.add_argument("--variants", required=True)
    ap.add_argument("--vcf", required=True, help="lifted VCF carrying INFO/R2 for those variants")
    ap.add_argument("--n-samples", type=int, required=True)
    ap.add_argument("--out", default=os.environ.get("EXP_OUT"))
    args = ap.parse_args()
    if not args.out:
        sys.exit("no output path (EXP_OUT unset and --out missing)")

    ann = pd.read_parquet(args.annotations, columns=["id", "gene_id", "MAF", "CADD_PHRED"])
    var = pd.read_parquet(args.variants, columns=["id"] + KEYS)
    var["chrom"] = var["chrom"].astype(str).str.replace("^chr", "", regex=True).astype("int64")

    r2 = vcf_r2(args.vcf)
    joined = var.merge(r2, on=KEYS, how="left")
    ann = ann.merge(joined[["id", "R2"]], on="id", how="left")

    maf = pd.to_numeric(ann["MAF"], errors="coerce")
    cadd = pd.to_numeric(ann["CADD_PHRED"], errors="coerce")
    r2v = pd.to_numeric(ann["R2"], errors="coerce")
    elig = (maf > 0) & (maf < MAF_MAX) & (cadd > CADD_MIN) & ann["gene_id"].notna()

    two_n = 2 * args.n_samples
    d = ann.loc[elig, ["gene_id"]].copy()
    d["mac"] = two_n * maf[elig].to_numpy()
    d["mac_eff"] = d["mac"] * np.nan_to_num(r2v[elig].to_numpy(), nan=0.0)
    g = d.groupby("gene_id").agg(n_variants=("mac", "size"),
                                 cmac=("mac", "sum"),
                                 cmac_eff=("mac_eff", "sum"))

    lines = []
    def emit(k, v):
        lines.append(f"{k}\t{v}")

    emit("annotation_rows", len(ann))
    emit("r2_join_nonnull_rows", int(r2v.notna().sum()))
    emit("r2_join_null_rows", int(r2v.isna().sum()))
    emit("eligible_rows", int(elig.sum()))
    emit("genes_with_eligible_variant", len(g))
    emit("n_samples", args.n_samples)
    emit("two_n", two_n)

    for name, col in (("cmac", "cmac"), ("cmac_eff", "cmac_eff"), ("n_variants", "n_variants")):
        s = g[col]
        for q, lab in ((0.05, "p05"), (0.25, "p25"), (0.5, "median"), (0.75, "p75"), (0.95, "p95")):
            emit(f"{name}_{lab}", f"{float(s.quantile(q)):.4f}")
        emit(f"{name}_min", f"{float(s.min()):.4f}")
        emit(f"{name}_max", f"{float(s.max()):.4f}")
        emit(f"{name}_mean", f"{float(s.mean()):.4f}")

    for thr in (0, 5, 10, 20, 50):
        emit(f"genes_pass_cmac_ge_{thr}", int((g["cmac"] >= thr).sum()))
        emit(f"genes_pass_cmac_eff_ge_{thr}", int((g["cmac_eff"] >= thr).sum()))
        emit(f"genes_lost_by_switching_to_eff_at_{thr}",
             int(((g["cmac"] >= thr) & (g["cmac_eff"] < thr)).sum()))
    ratio = (g["cmac_eff"] / g["cmac"].replace(0, np.nan)).dropna()
    for q, lab in ((0.05, "p05"), (0.5, "median"), (0.95, "p95")):
        emit(f"cmac_eff_over_cmac_{lab}", f"{float(ratio.quantile(q)):.4f}")
    emit("bonferroni_alpha_if_universe_equals_genes_with_eligible",
         f"{0.05/max(1,len(g)):.6g}")

    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write("key\tvalue\n")
        fh.write("\n".join(lines) + "\n")
    print(f"F3A_OK rows={len(lines)+1} genes={len(g)} out={args.out}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
