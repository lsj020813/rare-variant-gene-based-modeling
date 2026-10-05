#!/usr/bin/env python

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

CHI2_MED = stats.chi2.ppf(0.5, df=1)

def lam(p: np.ndarray) -> float:
    return float(np.median(stats.chi2.isf(p, df=1)) / CHI2_MED)

def ci(v: np.ndarray, lo=2.5, hi=97.5) -> list[float]:
    return [float(np.percentile(v, lo)), float(np.percentile(v, hi))]

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--new-results", required=True)
    ap.add_argument("--baseline-results", required=True)
    ap.add_argument("--saige-set", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--n-boot", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=20260728)
    ap.add_argument("--out-json", required=True)
    args = ap.parse_args()

    g = pd.read_parquet(args.genes, columns=["id", "gene", "gene_name"])
    ens = dict(zip(g["id"].astype("int64"), g["gene"]))

    n = pd.read_parquet(args.new_results)
    b = pd.read_parquet(args.baseline_results)
    s = pd.read_csv(args.saige_set, sep="\t").drop_duplicates(subset=["Region"])

    bmap = dict(zip(b["gene"].astype("int64"), b["pval"].astype(float)))
    smap = {r.Region: r for r in s.itertuples()}

    rows = []
    for r in n.itertuples():
        gid = int(r.gene)
        e = ens.get(gid)
        if gid in bmap and e in smap:
            sr = smap[e]
            rows.append((
                float(r.pval), bmap[gid],
                float(sr.Pvalue), float(sr.Pvalue_Burden), float(sr.Pvalue_SKAT),
            ))
    if not rows:
        raise SystemExit("JOIN FAILURE: no gene present in all three result sets")
    a = np.asarray(rows, dtype=float)
    p_gt, p_gp, sko, sbu, ska = a[:, 0], a[:, 1], a[:, 2], a[:, 3], a[:, 4]
    m = len(a)

    out: dict = {"genes_in_all_three": m, "n_boot": args.n_boot, "seed": args.seed}
    out["point"] = {
        "lambda_gc_gt": lam(p_gt),
        "lambda_gc_gp90": lam(p_gp),
        "delta_lambda_gc": lam(p_gt) - lam(p_gp),
    }

    nl_gt, nl_gp = -np.log10(p_gt), -np.log10(p_gp)
    sa = {"SKATO": -np.log10(sko), "Burden": -np.log10(sbu), "SKAT": -np.log10(ska)}
    for k, v in sa.items():
        out["point"][f"rho_after_{k}"] = float(stats.spearmanr(nl_gt, v).statistic)
        out["point"][f"rho_before_{k}"] = float(stats.spearmanr(nl_gp, v).statistic)
        out["point"][f"delta_rho_{k}"] = (
            out["point"][f"rho_after_{k}"] - out["point"][f"rho_before_{k}"]
        )

    rng = np.random.default_rng(args.seed)
    acc: dict[str, list[float]] = {k: [] for k in
                                   ["lam_gt", "lam_gp", "d_lam"]
                                   + [f"d_rho_{k}" for k in sa]
                                   + [f"rho_after_{k}" for k in sa]
                                   + [f"rho_before_{k}" for k in sa]}
    for _ in range(args.n_boot):
        idx = rng.integers(0, m, m)
        lg, lp = lam(p_gt[idx]), lam(p_gp[idx])
        acc["lam_gt"].append(lg)
        acc["lam_gp"].append(lp)
        acc["d_lam"].append(lg - lp)
        ngt, ngp = nl_gt[idx], nl_gp[idx]
        for k, v in sa.items():
            ra = float(stats.spearmanr(ngt, v[idx]).statistic)
            rb = float(stats.spearmanr(ngp, v[idx]).statistic)
            acc[f"rho_after_{k}"].append(ra)
            acc[f"rho_before_{k}"].append(rb)
            acc[f"d_rho_{k}"].append(ra - rb)

    out["bootstrap_ci95"] = {k: ci(np.asarray(v)) for k, v in acc.items()}
    out["bootstrap_p_two_sided_delta_ge_0"] = {
        k: float(2 * min((np.asarray(v) <= 0).mean(), (np.asarray(v) >= 0).mean()))
        for k, v in acc.items() if k.startswith("d_")
    }
    out["bootstrap_fraction_delta_positive"] = {
        k: float((np.asarray(v) > 0).mean()) for k, v in acc.items() if k.startswith("d_")
    }

    Path(args.out_json).write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
