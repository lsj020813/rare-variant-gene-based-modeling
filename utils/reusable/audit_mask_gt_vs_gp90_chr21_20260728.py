#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

MAF_CUT = 1e-3
CADD_CUT = 5.0

def mask_frame(anno: pd.DataFrame, maf_col: str) -> pd.DataFrame:
    return anno[(anno[maf_col] < MAF_CUT) & (anno["CADD_PHRED"] > CADD_CUT) & anno["gene_id"].notna()]

def describe(series: pd.Series) -> dict:
    s = pd.to_numeric(series, errors="coerce").dropna()
    if s.empty:
        return {"n": 0}
    qs = [0, 1, 5, 10, 25, 50, 75, 90, 95, 99, 100]
    out = {"n": int(len(s)), "mean": float(s.mean())}
    for q in qs:
        out[f"p{q}"] = float(np.percentile(s, q))
    return out

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotations-gp90", required=True, help="production annotation parquet (GP90 MAF)")
    ap.add_argument("--annotations-gt", required=True, help="new GT-frequency annotation parquet")
    ap.add_argument("--variants", required=True)
    ap.add_argument("--cohort-af-gp90", required=True, help="frozen GP90 cohort AF TSV (carries AC_GT/AN_GT)")
    ap.add_argument("--cohort-af-gt", required=True, help="new GT cohort AF TSV")
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--out-txt", required=True)
    args = ap.parse_args()

    out_json = Path(args.out_json)
    out_txt = Path(args.out_txt)
    for p in (out_json, out_txt):
        if p.exists():
            raise SystemExit(f"REFUSING_TO_OVERWRITE={p}")

    res: dict = {"thresholds": {"MAF": f"< {MAF_CUT}", "CADD_PHRED": f"> {CADD_CUT}"}}

    cols = ["id", "gene_id", "MAF", "MAF_MB", "AF", "CADD_PHRED"]
    a90 = pd.read_parquet(args.annotations_gp90, columns=cols)
    agt = pd.read_parquet(args.annotations_gt, columns=cols + ["AC_GT", "AN_GT"])

    res["annotation_rows_gp90"] = int(len(a90))
    res["annotation_rows_gt"] = int(len(agt))
    res["annotation_rows_identical"] = bool(len(a90) == len(agt))
    res["annotation_id_sequence_identical"] = bool(
        len(a90) == len(agt) and np.array_equal(a90["id"].to_numpy(), agt["id"].to_numpy())
    )
    res["cadd_column_identical"] = bool(
        res["annotation_id_sequence_identical"]
        and np.allclose(
            a90["CADD_PHRED"].to_numpy(dtype=float),
            agt["CADD_PHRED"].to_numpy(dtype=float),
            equal_nan=True,
        )
    )
    res["gene_id_column_identical"] = bool(
        res["annotation_id_sequence_identical"]
        and a90["gene_id"].fillna(-1).equals(agt["gene_id"].fillna(-1))
    )

    res["gp90_maf_constant_within_id"] = int(a90.groupby("id")["MAF"].nunique().max()) == 1
    res["gt_maf_constant_within_id"] = int(agt.groupby("id")["MAF"].nunique().max()) == 1

    key = ["chrom", "pos", "ref", "alt"]
    f90 = pd.read_csv(args.cohort_af_gp90, sep="\t", na_values=["NA"])
    fgt = pd.read_csv(args.cohort_af_gt, sep="\t", na_values=["NA"])
    for f in (f90, fgt):
        f["chrom"] = f["chrom"].astype(str).str.replace("^chr", "", regex=True)
        f["pos"] = f["pos"].astype("int64")
        f["ref"] = f["ref"].astype(str)
        f["alt"] = f["alt"].astype(str)
    left = f90[key + ["AC_GT", "AN_GT", "MAF_GT", "AC_GP90", "AN_GP90", "MAF"]].rename(
        columns={
            "AC_GT": "AC_GT_frozen",
            "AN_GT": "AN_GT_frozen",
            "MAF_GT": "MAF_GT_frozen",
            "MAF": "MAF_GP90",
        }
    )
    right = fgt[key + ["AC_GT", "AN_GT", "MAF"]].rename(
        columns={"AC_GT": "AC_GT_new", "AN_GT": "AN_GT_new", "MAF": "MAF_new"}
    )
    j = left.merge(right, on=key, how="inner", validate="1:1")
    res["crossval_rows_frozen"] = int(len(f90))
    res["crossval_rows_new"] = int(len(fgt))
    res["crossval_joined"] = int(len(j))
    res["crossval_AC_GT_identical"] = int((j["AC_GT_frozen"] == j["AC_GT_new"]).sum())
    res["crossval_AN_GT_identical"] = int((j["AN_GT_frozen"] == j["AN_GT_new"]).sum())
    res["crossval_MAF_max_abs_diff"] = float(
        np.abs(j["MAF_GT_frozen"].astype(float) - j["MAF_new"].astype(float)).max()
    )
    res["crossval_AC_GT_ge_AC_GP90_all"] = int((j["AC_GT_new"] >= j["AC_GP90"]).sum())

    m90 = mask_frame(a90, "MAF")
    mgt = mask_frame(agt, "MAF")

    def mask_stats(m: pd.DataFrame) -> dict:
        return {
            "annotation_rows_variant_gene_pairs": int(len(m)),
            "distinct_variant_gene_pairs": int(m[["id", "gene_id"]].drop_duplicates().shape[0]),
            "distinct_variants": int(m["id"].nunique()),
            "distinct_genes": int(m["gene_id"].nunique()),
            "variants_with_maf_exactly_zero": int(
                m.drop_duplicates("id").assign(z=lambda d: d["MAF"] == 0)["z"].sum()
            ),
            "maf_mb_rows_at_10000": int(np.isclose(m["MAF_MB"].to_numpy(dtype=float), 10000.0).sum()),
        }

    res["mask_gp90"] = mask_stats(m90)
    res["mask_gt"] = mask_stats(mgt)
    res["reference_expectation_gp90"] = {
        "distinct_variants": 4130,
        "annotation_rows_variant_gene_pairs": 4550,
        "distinct_genes": 212,
        "variants_with_maf_exactly_zero": 485,
    }
    res["gp90_mask_reproduces_reference"] = bool(
        res["mask_gp90"]["distinct_variants"] == 4130
        and res["mask_gp90"]["annotation_rows_variant_gene_pairs"] == 4550
        and res["mask_gp90"]["distinct_genes"] == 212
        and res["mask_gp90"]["variants_with_maf_exactly_zero"] == 485
    )

    ids90 = set(m90["id"])
    idsgt = set(mgt["id"])
    res["delta_variants"] = {
        "in_both": len(ids90 & idsgt),
        "left_gp90_only": len(ids90 - idsgt),
        "entered_gt_only": len(idsgt - ids90),
        "net_change": len(idsgt) - len(ids90),
    }
    g90 = set(m90["gene_id"].dropna())
    ggt = set(mgt["gene_id"].dropna())
    res["delta_genes"] = {
        "in_both": len(g90 & ggt),
        "left_gp90_only": len(g90 - ggt),
        "entered_gt_only": len(ggt - g90),
        "net_change": len(ggt) - len(g90),
    }
    p90 = set(map(tuple, m90[["id", "gene_id"]].drop_duplicates().to_numpy()))
    pgt = set(map(tuple, mgt[["id", "gene_id"]].drop_duplicates().to_numpy()))
    res["delta_variant_gene_pairs"] = {
        "in_both": len(p90 & pgt),
        "left_gp90_only": len(p90 - pgt),
        "entered_gt_only": len(pgt - p90),
        "net_change": len(pgt) - len(p90),
    }

    v90 = a90.drop_duplicates("id").set_index("id")[["MAF", "MAF_MB", "AF"]]
    vgt = agt.drop_duplicates("id").set_index("id")[["MAF", "MAF_MB", "AF", "AC_GT", "AN_GT"]]

    phantom = sorted(i for i in ids90 if v90.loc[i, "MAF"] == 0)
    res["gp90_phantom_n"] = len(phantom)
    ph_gt_maf = vgt.loc[phantom, "MAF"].astype(float)
    ph_gt_ac = vgt.loc[phantom, "AC_GT"].astype(int)
    res["phantoms"] = {
        "n": len(phantom),
        "gt_maf_gt0_no_longer_phantom": int((ph_gt_maf > 0).sum()),
        "gt_maf_still_exactly_zero": int((ph_gt_maf == 0).sum()),
        "gt_ac_gt0": int((ph_gt_ac > 0).sum()),
        "gt_ac_eq0": int((ph_gt_ac == 0).sum()),
        "still_in_gt_mask": len([i for i in phantom if i in idsgt]),
        "left_gt_mask": len([i for i in phantom if i not in idsgt]),
    }

    gp90_ids = sorted(ids90)
    gt_maf_of_gp90 = vgt.loc[gp90_ids, "MAF"].astype(float)
    would_exclude = [i for i, m in gt_maf_of_gp90.items() if not (m < MAF_CUT)]
    res["gt_cut_would_exclude"] = {
        "n": len(would_exclude),
        "reference_expectation": 321,
        "actually_left_gt_mask": len([i for i in would_exclude if i not in idsgt]),
        "still_in_gt_mask": len([i for i in would_exclude if i in idsgt]),
        "gt_maf_max_among_them": float(gt_maf_of_gp90.loc[would_exclude].max()) if would_exclude else None,
    }
    left = sorted(ids90 - idsgt)
    left_gt_maf = vgt.loc[left, "MAF"].astype(float) if left else pd.Series(dtype=float)
    res["left_decomposition"] = {
        "n_left": len(left),
        "left_because_gt_maf_ge_cut": int((left_gt_maf >= MAF_CUT).sum()) if len(left) else 0,
        "left_for_other_reasons": int((left_gt_maf < MAF_CUT).sum()) if len(left) else 0,
    }

    entered = sorted(idsgt - ids90)
    res["entrants"] = {"n": len(entered)}
    if entered:
        e90 = v90.loc[entered]
        egt = vgt.loc[entered]
        res["entrants"].update(
            {
                "gp90_maf_ge_cut": int((e90["MAF"].astype(float) >= MAF_CUT).sum()),
                "gp90_maf_lt_cut_but_excluded_for_other_reason": int(
                    (e90["MAF"].astype(float) < MAF_CUT).sum()
                ),
                "gt_maf_quantiles": describe(egt["MAF"]),
                "gp90_maf_quantiles": describe(e90["MAF"]),
                "gt_ac_eq0": int((egt["AC_GT"].astype(int) == 0).sum()),
                "distinct_genes_touched": int(
                    mgt[mgt["id"].isin(entered)]["gene_id"].nunique()
                ),
                "genes_new_to_the_mask_because_of_them": len(
                    set(mgt[mgt["id"].isin(entered)]["gene_id"].dropna()) - g90
                ),
            }
        )

    res["maf_mb"] = {
        "mask_gp90_variant_level": describe(v90.loc[gp90_ids, "MAF_MB"]),
        "mask_gt_variant_level": describe(vgt.loc[sorted(idsgt), "MAF_MB"]),
        "mask_gp90_variant_level_at_10000": int(
            np.isclose(v90.loc[gp90_ids, "MAF_MB"].to_numpy(dtype=float), 10000.0).sum()
        ),
        "mask_gt_variant_level_at_10000": int(
            np.isclose(vgt.loc[sorted(idsgt), "MAF_MB"].to_numpy(dtype=float), 10000.0).sum()
        ),
        "mask_gp90_annotation_rows_at_10000": res["mask_gp90"]["maf_mb_rows_at_10000"],
        "mask_gt_annotation_rows_at_10000": res["mask_gt"]["maf_mb_rows_at_10000"],
        "chromwide_gp90_annotation_rows_at_10000": int(
            np.isclose(a90["MAF_MB"].to_numpy(dtype=float), 10000.0).sum()
        ),
        "chromwide_gt_annotation_rows_at_10000": int(
            np.isclose(agt["MAF_MB"].to_numpy(dtype=float), 10000.0).sum()
        ),
        "note": (
            "10000 == (1e-8)**-0.5 exactly, so it is a COMPUTED value for AF==0 and is "
            "numerically indistinguishable from the upstream missing-value fill constant."
        ),
    }
    both = sorted(ids90 & idsgt)
    if both:
        d = vgt.loc[both, "MAF_MB"].astype(float) - v90.loc[both, "MAF_MB"].astype(float)
        res["maf_mb"]["paired_delta_on_common_variants"] = describe(d)
        res["maf_mb"]["paired_ratio_on_common_variants"] = describe(
            vgt.loc[both, "MAF_MB"].astype(float) / v90.loc[both, "MAF_MB"].astype(float)
        )

    res["chromwide"] = {
        "variants": int(len(j)),
        "maf_ratio_gt_over_gp90_quantiles": describe(
            (j["MAF_new"].astype(float) / j["MAF_GP90"].astype(float)).replace([np.inf, -np.inf], np.nan)
        ),
        "an_gp90_lt_an_gt": int((j["AN_GP90"] < j["AN_GT_new"]).sum()),
        "maf_gt_ge_maf_gp90": int((j["MAF_new"].astype(float) >= j["MAF_GP90"].astype(float)).sum()),
        "ac_gt_eq0": int((j["AC_GT_new"] == 0).sum()),
        "ac_gp90_eq0": int((j["AC_GP90"] == 0).sum()),
    }

    def gene_variant_counts(m: pd.DataFrame) -> pd.Series:
        return m.drop_duplicates(["id", "gene_id"]).groupby("gene_id")["id"].size()

    c90 = gene_variant_counts(m90)
    cgt = gene_variant_counts(mgt)
    res["per_gene_variant_count"] = {
        "gp90": describe(c90),
        "gt": describe(cgt),
        "genes_losing_all_variants": int(len(set(c90.index) - set(cgt.index))),
        "genes_gaining_entry": int(len(set(cgt.index) - set(c90.index))),
    }

    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(res, indent=2, sort_keys=True))

    lines = [
        "MASK DELTA DIGEST  GT-derived cohort frequency vs GP>=0.90  chr21  20260728",
        f"thresholds unchanged: MAF < {MAF_CUT}, CADD_PHRED > {CADD_CUT}",
        "",
        f"GP90 mask reproduces the frozen reference: {res['gp90_mask_reproduces_reference']}",
        f"  GP90  variants={res['mask_gp90']['distinct_variants']} "
        f"pairs(annotation rows)={res['mask_gp90']['annotation_rows_variant_gene_pairs']} "
        f"pairs(distinct)={res['mask_gp90']['distinct_variant_gene_pairs']} "
        f"genes={res['mask_gp90']['distinct_genes']}",
        f"  GT    variants={res['mask_gt']['distinct_variants']} "
        f"pairs(annotation rows)={res['mask_gt']['annotation_rows_variant_gene_pairs']} "
        f"pairs(distinct)={res['mask_gt']['distinct_variant_gene_pairs']} "
        f"genes={res['mask_gt']['distinct_genes']}",
        "",
        f"variants: both={res['delta_variants']['in_both']} left={res['delta_variants']['left_gp90_only']} "
        f"entered={res['delta_variants']['entered_gt_only']} net={res['delta_variants']['net_change']}",
        f"genes:    both={res['delta_genes']['in_both']} left={res['delta_genes']['left_gp90_only']} "
        f"entered={res['delta_genes']['entered_gt_only']} net={res['delta_genes']['net_change']}",
        "",
        f"phantoms (GP90 MAF==0): n={res['phantoms']['n']} "
        f"no_longer_phantom={res['phantoms']['gt_maf_gt0_no_longer_phantom']} "
        f"still_zero={res['phantoms']['gt_maf_still_exactly_zero']} "
        f"left_mask={res['phantoms']['left_gt_mask']}",
        f"GT cut would exclude: n={res['gt_cut_would_exclude']['n']} "
        f"(reference 321) actually_left={res['gt_cut_would_exclude']['actually_left_gt_mask']}",
        f"entrants: n={res['entrants']['n']}",
        "",
        f"MAF_MB at 10000 (variant level): GP90={res['maf_mb']['mask_gp90_variant_level_at_10000']} "
        f"GT={res['maf_mb']['mask_gt_variant_level_at_10000']}",
        f"MAF_MB at 10000 (annotation rows): GP90={res['maf_mb']['mask_gp90_annotation_rows_at_10000']} "
        f"GT={res['maf_mb']['mask_gt_annotation_rows_at_10000']}",
        "",
        f"crossval vs frozen table: joined={res['crossval_joined']} "
        f"AC_GT identical={res['crossval_AC_GT_identical']} AN_GT identical={res['crossval_AN_GT_identical']} "
        f"MAF max abs diff={res['crossval_MAF_max_abs_diff']}",
    ]
    out_txt.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"WROTE {out_json}")
    print(f"WROTE {out_txt}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
