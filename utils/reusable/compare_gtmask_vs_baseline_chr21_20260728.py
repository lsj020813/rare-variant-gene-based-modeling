#!/usr/bin/env python

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

CHI2_1DF_MEDIAN = stats.chi2.ppf(0.5, df=1)

def lambda_gc(p: np.ndarray) -> float:
    p = np.asarray(p, dtype=float)
    p = p[np.isfinite(p) & (p > 0)]
    if not len(p):
        return float("nan")
    return float(np.median(stats.chi2.isf(p, df=1)) / CHI2_1DF_MEDIAN)

def bh_qvalues(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.minimum(q, 1.0)
    return out

def arm_summary(p: np.ndarray, m_genome: int = 20000) -> dict:
    p = np.asarray(p, dtype=float)
    p = p[np.isfinite(p)]
    m = len(p)
    q = bh_qvalues(p) if m else np.array([])
    return {
        "genes_tested": int(m),
        "min_p": float(p.min()) if m else None,
        "median_p": float(np.median(p)) if m else None,
        "lambda_gc": lambda_gc(p),
        "n_nominal_0.05": int((p < 0.05).sum()),
        "expected_nominal_0.05": float(0.05 * m),
        "bonferroni_threshold_over_tested": float(0.05 / m) if m else None,
        "n_bonferroni_over_tested": int((p < 0.05 / m).sum()) if m else 0,
        "n_bonferroni_over_20000": int((p < 0.05 / m_genome).sum()),
        "n_bh_0.05": int((q <= 0.05).sum()) if m else 0,
        "n_bh_0.10": int((q <= 0.10).sum()) if m else 0,
        "min_bh_q": float(q.min()) if m else None,
        "ks_vs_uniform_D": float(stats.kstest(p, "uniform").statistic) if m else None,
        "ks_vs_uniform_p": float(stats.kstest(p, "uniform").pvalue) if m else None,
    }

def load_results(path: str, genes: pd.DataFrame, label: str, out: dict) -> pd.DataFrame:
    d = pd.read_parquet(path)
    out[f"{label}_rows"] = int(len(d))
    out[f"{label}_columns"] = list(d.columns)
    out[f"{label}_methods"] = sorted(d["Method"].astype(str).unique().tolist())
    out[f"{label}_dup_gene"] = int(d.duplicated(["gene"]).sum())
    out[f"{label}_dup_gene_method"] = int(d.duplicated(["gene", "Method"]).sum())
    out[f"{label}_dtype_gene_raw"] = str(d["gene"].dtype)
    d = d[["gene", "pval", "beta"]].copy()
    d["gene"] = d["gene"].astype("int64")
    assert str(d["gene"].dtype) == "int64"
    m = d.merge(genes, left_on="gene", right_on="id", how="left", validate="1:1")
    out[f"{label}_bridged_to_ensembl"] = int(m["ensembl"].notna().sum())
    out[f"{label}_bridge_hit_rate"] = float(m["ensembl"].notna().mean())
    if m["ensembl"].notna().sum() == 0:
        raise SystemExit(f"JOIN FAILURE: 0 of {len(m)} {label} genes bridged to Ensembl")
    m["gene"] = m["gene"].astype("int64")
    return m.rename(columns={"pval": f"pval_{label}", "beta": f"beta_{label}"})

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--new-results", required=True)
    ap.add_argument("--baseline-results", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--out-txt", required=True)
    args = ap.parse_args()

    for p in (Path(args.out_json), Path(args.out_txt)):
        if p.exists():
            raise SystemExit(f"REFUSING_TO_OVERWRITE={p}")

    out: dict = {
        "new_results": args.new_results,
        "baseline_results": args.baseline_results,
    }

    g = pd.read_parquet(args.genes, columns=["id", "gene", "gene_name"])
    out["gene_table_rows"] = int(len(g))
    out["gene_table_dup_id"] = int(g.duplicated(["id"]).sum())
    out["gene_table_dup_ensembl"] = int(g.duplicated(["gene"]).sum())
    assert out["gene_table_dup_id"] == 0 and out["gene_table_dup_ensembl"] == 0
    g = g.rename(columns={"gene": "ensembl"})
    g["id"] = g["id"].astype("int64")

    new = load_results(args.new_results, g, "gt", out)
    base = load_results(args.baseline_results, g, "gp90", out)

    j_int = new.merge(base[["gene", "pval_gp90", "beta_gp90"]], on="gene",
                      how="inner", validate="1:1")
    j_ens = new.merge(base[["ensembl", "pval_gp90", "beta_gp90"]], on="ensembl",
                      how="inner", validate="1:1")
    out["paired_on_integer_gene_id"] = int(len(j_int))
    out["paired_on_ensembl_id"] = int(len(j_ens))
    out["two_pairings_identical"] = bool(
        len(j_int) == len(j_ens)
        and set(j_int["ensembl"]) == set(j_ens["ensembl"])
    )
    if len(j_int) == 0:
        raise SystemExit("JOIN FAILURE: 0 genes shared between the two runs")
    out["join_hit_rate_new_side"] = float(len(j_int) / len(new))
    out["join_hit_rate_baseline_side"] = float(len(j_int) / len(base))

    only_new = sorted(set(new["ensembl"].dropna()) - set(base["ensembl"].dropna()))
    only_base = sorted(set(base["ensembl"].dropna()) - set(new["ensembl"].dropna()))
    sym = dict(zip(g["ensembl"], g["gene_name"]))
    out["genes_only_in_gt_run"] = {
        "n": len(only_new),
        "symbols": sorted(str(sym.get(e)) for e in only_new),
    }
    out["genes_only_in_gp90_baseline"] = {
        "n": len(only_base),
        "symbols": sorted(str(sym.get(e)) for e in only_base),
    }

    j = j_int
    out["arm_gt"] = arm_summary(new["pval_gt"].to_numpy())
    out["arm_gp90"] = arm_summary(base["pval_gp90"].to_numpy())
    out["arm_gt_on_paired"] = arm_summary(j["pval_gt"].to_numpy())
    out["arm_gp90_on_paired"] = arm_summary(j["pval_gp90"].to_numpy())
    out["delta_lambda_gc_full"] = (
        out["arm_gt"]["lambda_gc"] - out["arm_gp90"]["lambda_gc"]
    )
    out["delta_lambda_gc_paired"] = (
        out["arm_gt_on_paired"]["lambda_gc"] - out["arm_gp90_on_paired"]["lambda_gc"]
    )

    pn = j["pval_gt"].to_numpy(dtype=float)
    pb = j["pval_gp90"].to_numpy(dtype=float)
    nl_n, nl_b = -np.log10(pn), -np.log10(pb)
    ok = np.isfinite(nl_n) & np.isfinite(nl_b)
    sp = stats.spearmanr(nl_n[ok], nl_b[ok])
    kd = stats.kendalltau(nl_n[ok], nl_b[ok])
    pe = stats.pearsonr(nl_n[ok], nl_b[ok])
    out["neglog10p_agreement"] = {
        "n": int(ok.sum()),
        "spearman_rho": float(sp.statistic),
        "spearman_p": float(sp.pvalue),
        "kendall_tau": float(kd.statistic),
        "pearson_r": float(pe.statistic),
    }

    dlog = np.abs(nl_n - nl_b)
    dp = np.abs(pn - pb)
    out["delta"] = {
        "max_abs_delta_log10p": float(np.nanmax(dlog)),
        "median_abs_delta_log10p": float(np.nanmedian(dlog)),
        "p90_abs_delta_log10p": float(np.nanpercentile(dlog, 90)),
        "n_genes_abs_delta_p_gt_0.01": int((dp > 0.01).sum()),
        "n_genes_abs_delta_p_gt_0.05": int((dp > 0.05).sum()),
        "n_genes_p_identical": int((pn == pb).sum()),
        "max_abs_delta_p": float(np.nanmax(dp)),
        "median_abs_delta_p": float(np.nanmedian(dp)),
    }
    imax = int(np.nanargmax(dlog))
    out["delta"]["largest_mover"] = {
        "gene_name": str(j["gene_name"].iloc[imax]),
        "p_gp90": float(pb[imax]),
        "p_gt": float(pn[imax]),
        "abs_delta_log10p": float(dlog[imax]),
    }

    bn = j["beta_gt"].to_numpy(dtype=float)
    bb = j["beta_gp90"].to_numpy(dtype=float)
    flip = np.isfinite(bn) & np.isfinite(bb) & (np.sign(bn) != np.sign(bb))
    out["beta"] = {
        "n_compared": int((np.isfinite(bn) & np.isfinite(bb)).sum()),
        "n_sign_flips": int(flip.sum()),
        "sign_flip_symbols": sorted(j["gene_name"][flip].astype(str).tolist()),
        "spearman_rho": float(stats.spearmanr(bn, bb).statistic),
        "pearson_r": float(stats.pearsonr(bn, bb).statistic),
        "max_abs_delta_beta": float(np.nanmax(np.abs(bn - bb))),
        "median_abs_delta_beta": float(np.nanmedian(np.abs(bn - bb))),
        "n_beta_positive_gt": int((bn > 0).sum()),
        "n_beta_positive_gp90": int((bb > 0).sum()),
    }

    for symbol in ("CBS", "SLC19A1"):
        for tag, frame, pcol in (("gt", new, "pval_gt"), ("gp90", base, "pval_gp90")):
            f = frame.dropna(subset=[pcol]).copy()
            f["_r"] = f[pcol].rank(method="min")
            row = f[f["gene_name"] == symbol]
            out.setdefault("positive_controls", {}).setdefault(symbol, {})[tag] = (
                None if row.empty else {
                    "rank": int(row["_r"].iloc[0]),
                    "of": int(len(f)),
                    "p": float(row[pcol].iloc[0]),
                    "beta": float(row[f"beta_{tag}"].iloc[0]),
                }
            )

    top_gt = new.nsmallest(10, "pval_gt")[["gene_name", "pval_gt", "beta_gt"]]
    out["gt_top10"] = [
        {"gene_name": str(r.gene_name), "p": float(r.pval_gt), "beta": float(r.beta_gt)}
        for r in top_gt.itertuples()
    ]
    top_bp = base.nsmallest(10, "pval_gp90")[["gene_name", "pval_gp90"]]
    out["gp90_top10"] = [
        {"gene_name": str(r.gene_name), "p": float(r.pval_gp90)}
        for r in top_bp.itertuples()
    ]

    Path(args.out_json).write_text(json.dumps(out, indent=2, default=str) + "\n")
    lines = [
        "chr21 HOMOCYST_rint DeepRVAT: GT-frequency mask vs GP90-frequency baseline",
        f"new       : {args.new_results}",
        f"baseline  : {args.baseline_results}",
        "",
        f"genes tested   GT {out['arm_gt']['genes_tested']}   GP90 {out['arm_gp90']['genes_tested']}",
        f"paired genes   {out['paired_on_integer_gene_id']} (int key) / "
        f"{out['paired_on_ensembl_id']} (Ensembl key), identical={out['two_pairings_identical']}",
        f"only in GT run     {out['genes_only_in_gt_run']['n']} {out['genes_only_in_gt_run']['symbols']}",
        f"only in baseline   {out['genes_only_in_gp90_baseline']['n']} {out['genes_only_in_gp90_baseline']['symbols']}",
        "",
        f"lambda_GC   GT {out['arm_gt']['lambda_gc']:.4f}   GP90 {out['arm_gp90']['lambda_gc']:.4f}   "
        f"delta {out['delta_lambda_gc_full']:+.4f}",
        f"min p       GT {out['arm_gt']['min_p']:.6g}   GP90 {out['arm_gp90']['min_p']:.6g}",
        f"Bonferroni over tested   GT {out['arm_gt']['n_bonferroni_over_tested']}   "
        f"GP90 {out['arm_gp90']['n_bonferroni_over_tested']}",
        f"BH 0.05 / 0.10   GT {out['arm_gt']['n_bh_0.05']}/{out['arm_gt']['n_bh_0.10']}   "
        f"GP90 {out['arm_gp90']['n_bh_0.05']}/{out['arm_gp90']['n_bh_0.10']}",
        "",
        f"Spearman rho of -log10 p across paired genes  {out['neglog10p_agreement']['spearman_rho']:+.4f} "
        f"(p={out['neglog10p_agreement']['spearman_p']:.3g})",
        f"max |delta log10 p|   {out['delta']['max_abs_delta_log10p']:.4f}",
        f"median |delta log10 p| {out['delta']['median_abs_delta_log10p']:.4f}",
        f"genes with |delta p| > 0.01   {out['delta']['n_genes_abs_delta_p_gt_0.01']}",
        f"genes with identical p        {out['delta']['n_genes_p_identical']}",
        f"beta sign flips   {out['beta']['n_sign_flips']} of {out['beta']['n_compared']}",
    ]
    Path(args.out_txt).write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
