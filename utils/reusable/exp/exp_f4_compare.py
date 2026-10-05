#!/usr/bin/env python
import sys

import numpy as np
import pandas as pd

VERDICT_EQUIVALENT = "PASS_PREFILTER_EQUIVALENT"
VERDICT_DIFFERS = "DIFFERS_WITHIN_PIPELINE_PASS"

def main(argv):
    mask_p, full_p, out_p, tag = argv[1:5]
    tol_p, tol_logp = float(argv[5]), float(argv[6])
    report_tols = [float(x) for x in argv[7].split()]

    lines = [
        f"tag\t{tag}",
        f"declared_tol_p\t{tol_p:g}",
        f"declared_tol_neglog10p\t{tol_logp:g}",
        "tolerance_units\tp_scale_and_neglog10p_scale_only; annotation-feature tolerances are not convertible",
        "equivalence_contract\tall_four_of(genes_exceeding_declared_tol_p,genes_exceeding_declared_tol_neglog10p,"
        "genes_only_in_full,genes_only_in_mask)_must_be_zero",
    ]

    def emit_fail(detail, extra=()):
        lines.extend([
            "pipeline_validation\tFAIL",
            f"pipeline_validation_detail\t{detail}",
            "verdict\tNOT_EVALUATED",
            "verdict_reason\tpipeline_validation_failed_no_scientific_claim_made",
        ])
        lines.extend(extra)
        with open(out_p, "w") as fh:
            fh.write("key\tvalue\n" + "\n".join(lines) + "\n")
        return 1

    try:
        m = pd.read_parquet(mask_p)
        f = pd.read_parquet(full_p)
    except Exception as e:
        return emit_fail(f"unreadable_input:{type(e).__name__}", ("genes_compared\t0",))

    problems = [
        f"{name}_missing_column_{col}"
        for name, d in (("mask", m), ("full", f))
        for col in ("gene", "pval")
        if col not in d.columns
    ]
    if problems:
        return emit_fail(";".join(problems), ("genes_compared\t0",))

    m["gene"] = m["gene"].astype("int64")
    f["gene"] = f["gene"].astype("int64")
    dup_m = int(m["gene"].duplicated().sum())
    dup_f = int(f["gene"].duplicated().sum())
    cols = [c for c in ("gene", "pval", "beta") if c in m.columns and c in f.columns]
    M = m[cols].drop_duplicates("gene").set_index("gene")
    F = f[cols].drop_duplicates("gene").set_index("gene")
    only_full = sorted(F.index.difference(M.index))
    only_mask = sorted(M.index.difference(F.index))
    n_only_full, n_only_mask = len(only_full), len(only_mask)
    universe_identical = (n_only_full == 0 and n_only_mask == 0)
    J = M.join(F, how="inner", lsuffix="_mask", rsuffix="_full")

    counts = [
        f"genes_full\t{len(F)}",
        f"genes_mask\t{len(M)}",
        f"genes_compared\t{len(J)}",
        f"genes_only_in_full\t{n_only_full}",
        f"genes_only_in_mask\t{n_only_mask}",
        f"dup_genes_full\t{dup_f}",
        f"dup_genes_mask\t{dup_m}",
        f"gene_universe_identical\t{universe_identical}",
        "key_universe\tp_values_are_compared_on_the_intersection; the gene universes themselves are "
        "compared separately and a difference blocks the equivalence verdict",
    ]
    lines.extend(counts)

    problems = []
    if len(J) == 0:
        problems.append("no_overlapping_genes")
    if dup_m or dup_f:
        problems.append(f"duplicate_gene_keys(mask={dup_m},full={dup_f})")
    pm = pf = None
    if len(J):
        pm = pd.to_numeric(J["pval_mask"], errors="coerce")
        pf = pd.to_numeric(J["pval_full"], errors="coerce")
        if not np.isfinite(pm).all() or not np.isfinite(pf).all():
            problems.append("non_finite_pvalues")
        elif ((pm <= 0) | (pm > 1) | (pf <= 0) | (pf > 1)).any():
            problems.append("pvalues_outside_(0,1]")

    if problems:
        na = [f"genes_exceeding_delta_neglog10p_{t:g}\tNA" for t in report_tols]
        na += [
            "genes_exceeding_declared_tol_p\tNA",
            "genes_exceeding_declared_tol_neglog10p\tNA",
            "max_abs_delta_pval\tNA",
            "max_abs_delta_neglog10p\tNA",
        ]
        return emit_fail(";".join(problems), na)

    detail = ["no_duplicate_gene_keys", "all_intersection_pvalues_finite_in_(0,1]"]
    detail.append(
        "identical_key_universe" if universe_identical
        else f"key_universe_DIFFERS(only_in_full={n_only_full},only_in_mask={n_only_mask})"
    )
    lines.append("pipeline_validation\tPASS")
    lines.append("pipeline_validation_detail\t" + ";".join(detail))

    dp = np.abs(pm - pf)
    lm, lf = -np.log10(pm), -np.log10(pf)
    dl = np.abs(lm - lf)
    n_p = int((dp > tol_p).sum())
    n_l = int((dl > tol_logp).sum())
    lines += [
        f"max_abs_delta_pval\t{float(dp.max()):.6e}",
        f"max_abs_delta_neglog10p\t{float(dl.max()):.6e}",
        f"median_abs_delta_neglog10p\t{float(np.median(dl)):.6e}",
        f"genes_exceeding_declared_tol_p\t{n_p}",
        f"genes_exceeding_declared_tol_neglog10p\t{n_l}",
    ]
    for t in report_tols:
        lines.append(f"genes_exceeding_delta_neglog10p_{t:g}\t{int((dl > t).sum())}")
    if len(J) > 2:
        lines.append(f"pearson_neglog10p\t{float(np.corrcoef(lm, lf)[0, 1]):.6f}")
    if "beta_mask" in J.columns and "beta_full" in J.columns:
        bm = pd.to_numeric(J["beta_mask"], errors="coerce")
        bf = pd.to_numeric(J["beta_full"], errors="coerce")
        lines.append(f"sign_flips_beta\t{int((np.sign(bm) != np.sign(bf)).sum())}")

    reasons = []
    if n_p:
        reasons.append(f"genes_exceeding_declared_tol_p={n_p}")
    if n_l:
        reasons.append(f"genes_exceeding_declared_tol_neglog10p={n_l}")
    if n_only_full:
        reasons.append(f"genes_only_in_full={n_only_full}")
    if n_only_mask:
        reasons.append(f"genes_only_in_mask={n_only_mask}")
    equivalent = not reasons
    lines.append("verdict\t" + (VERDICT_EQUIVALENT if equivalent else VERDICT_DIFFERS))
    lines.append("verdict_reason\t" + ("all_four_equivalence_counts_zero" if equivalent else ";".join(reasons)))
    lines.append(
        "caveat\tthe two arms also differ by annotation batch partition; F-2 quantifies that "
        "component in feature units, which is why it is not used as a p-value tolerance"
    )
    with open(out_p, "w") as fh:
        fh.write("key\tvalue\n" + "\n".join(lines) + "\n")
    print("F4_COMPARE", [l for l in lines if l.startswith(("verdict\t", "pipeline_validation\t"))])
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
