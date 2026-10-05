#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import os
from collections import defaultdict

import numpy as np
import pandas as pd
import yaml

KEYS = ["chrom", "pos", "ref", "alt"]

def load_arm(ann_path: str, var_path: str, feats: list[str]) -> tuple[pd.DataFrame, int]:
    ann = pd.read_parquet(ann_path)
    var = pd.read_parquet(var_path, columns=["id"] + KEYS)
    var["chrom"] = var["chrom"].astype(str).str.replace("^chr", "", regex=True)
    m = ann.merge(var, on="id", how="left", validate="m:1", sort=False)
    unresolved = int(m["pos"].isna().sum())
    m = m.loc[m["pos"].notna()].copy()
    m["_vh"] = [hashlib.sha1(f"{c}:{int(p)}:{r}:{a}".encode()).hexdigest()[:16]
                for c, p, r, a in zip(m["chrom"], m["pos"], m["ref"], m["alt"])]
    m["_gene"] = m["gene_id"].astype("Int64").astype(str)
    keep = ["_vh", "_gene", "gene_id"] + [f for f in feats if f in m.columns]
    return m[keep], unresolved

def numeric_matrix(frame: pd.DataFrame, feats: list[str]) -> tuple[np.ndarray, int]:
    raw = frame[feats]
    conv = raw.apply(pd.to_numeric, errors="coerce")
    coerced = int((conv.isna() & raw.notna()).to_numpy().sum())
    return conv.to_numpy(dtype=float), coerced

def canon_key(row: np.ndarray) -> tuple:
    return tuple((1, 0.0) if np.isnan(v) else (0, 0.0 if v == 0 else float(v)) for v in row)

def canonicalize(frame: pd.DataFrame, feats: list[str], group_col: str) -> tuple[pd.DataFrame, int]:
    mat, coerced = numeric_matrix(frame, feats)
    grouped: dict[str, list[np.ndarray]] = defaultdict(list)
    for g, row in zip(frame[group_col].astype(str), mat):
        grouped[g].append(row)
    idx: list[str] = []
    rows: list[np.ndarray] = []
    for g in sorted(grouped):
        for occ, row in enumerate(sorted(grouped[g], key=canon_key)):
            idx.append(f"{g}#{occ}")
            rows.append(row)
    out = pd.DataFrame(rows, columns=feats, index=pd.Index(idx, name="canonical_row_key"))
    return out, coerced

def compare(W: pd.DataFrame, S: pd.DataFrame, feats: list[str], abs_tol: float, tag: str,
            emit) -> tuple[int, int, float, float]:
    diff_cols = 0
    cells = 0
    worst_abs = 0.0
    worst_rel = 0.0
    for c in feats:
        a = W[c].to_numpy(dtype=float)
        b = S[c].to_numpy(dtype=float)
        equal = np.isclose(a, b, rtol=0.0, atol=abs_tol, equal_nan=True)
        nz = np.flatnonzero(~equal)
        if nz.size:
            diff_cols += 1
            cells += int(nz.size)
            finite = np.isfinite(a[nz]) & np.isfinite(b[nz])
            d = np.full(nz.size, np.inf)
            d[finite] = np.abs(a[nz][finite] - b[nz][finite])
            sc = np.maximum(np.abs(a[nz]), np.abs(b[nz]))
            sc[sc == 0] = 1.0
            fin = np.isfinite(d)
            if fin.any():
                worst_abs = max(worst_abs, float(d[fin].max()))
                worst_rel = max(worst_rel, float((d[fin] / sc[fin]).max()))
            else:
                worst_abs = float("inf")
                worst_rel = float("inf")
            emit(f"{tag}_col_{c}_differing_rows", int(nz.size))
            emit(f"{tag}_col_{c}_max_abs", f"{worst_abs:.6e}")
    return diff_cols, cells, worst_abs, worst_rel

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--whole-ann", required=True); ap.add_argument("--whole-var", required=True)
    ap.add_argument("--a-ann", required=True);     ap.add_argument("--a-var", required=True)
    ap.add_argument("--b-ann", required=True);     ap.add_argument("--b-var", required=True)
    ap.add_argument("--model-config", required=True)
    ap.add_argument("--abs-tol", type=float, default=0.0)
    ap.add_argument("--expect-whole-variants", type=int, default=8000)
    ap.add_argument("--expect-half-variants", type=int, default=4000)
    ap.add_argument("--out", default=os.environ.get("EXP_OUT"), required=False)
    args = ap.parse_args()
    if not args.out:
        raise SystemExit("no output path")

    feats_all = list(yaml.safe_load(open(args.model_config))["rare_variant_annotations"])
    W, uw = load_arm(args.whole_ann, args.whole_var, feats_all)
    A, ua = load_arm(args.a_ann, args.a_var, feats_all)
    B, ub = load_arm(args.b_ann, args.b_var, feats_all)
    feats = [f for f in feats_all if f in W.columns and f in A.columns and f in B.columns]

    lines: list[str] = []
    def emit(k, v):
        lines.append(f"{k}\t{v}")

    emit("validator", os.path.basename(__file__))
    emit("abs_tol_declared", args.abs_tol)
    emit("comparison_contract", "two_layer__L1_gene_assigned_basekey_multiset__L2_gene_unassigned_per_variant_multiset")
    emit("cross_arm_key", "sha1(chrom:pos:ref:alt)[:16]; per-arm id is an arm-local index and is NOT comparable across arms")
    emit("row_aggregation_applied", "False")
    emit("features_declared", len(feats_all))
    emit("features_present_in_all_arms", len(feats))
    emit("whole_unresolved_ids", uw); emit("a_unresolved_ids", ua); emit("b_unresolved_ids", ub)
    for tag, D in (("whole", W), ("a", A), ("b", B)):
        emit(f"{tag}_rows", len(D))
        emit(f"{tag}_variants", int(D["_vh"].nunique()))
        emit(f"{tag}_rows_gene_assigned", int(D["gene_id"].notna().sum()))
        emit(f"{tag}_rows_gene_unassigned", int(D["gene_id"].isna().sum()))

    S = pd.concat([A, B], ignore_index=True)

    vw, va, vb = set(W["_vh"]), set(A["_vh"]), set(B["_vh"])
    emit("variant_union_equals_whole", str((va | vb) == vw))
    emit("variant_intersection_empty", str(len(va & vb) == 0))
    emit("variant_counts", f"whole={len(vw)} a={len(va)} b={len(vb)}")
    counts_ok = (len(vw) == args.expect_whole_variants and len(va) == args.expect_half_variants
                 and len(vb) == args.expect_half_variants)
    emit("expected_variant_counts", f"whole={args.expect_whole_variants} a=b={args.expect_half_variants}")
    emit("variant_counts_as_designed", str(counts_ok))

    fatal = (uw or ua or ub or (va | vb) != vw or len(va & vb) or not counts_ok)

    gw = W.loc[W["gene_id"].notna()].copy(); gs = S.loc[S["gene_id"].notna()].copy()
    gw["_k"] = gw["_vh"] + ":" + gw["_gene"]; gs["_k"] = gs["_vh"] + ":" + gs["_gene"]
    emit("l1_rows_whole", len(gw)); emit("l1_rows_split", len(gs))
    emit("l1_dup_keys_whole", int(gw["_k"].duplicated().sum()))
    emit("l1_dup_keys_split", int(gs["_k"].duplicated().sum()))
    kw, ks = set(gw["_k"]), set(gs["_k"])
    emit("l1_basekeys_whole", len(kw)); emit("l1_basekeys_split", len(ks))
    emit("l1_basekeys_only_in_whole", len(kw - ks))
    emit("l1_basekeys_only_in_split", len(ks - kw))
    wc = gw["_k"].value_counts().to_dict(); sc = gs["_k"].value_counts().to_dict()
    l1_mult_mismatch = sum(wc.get(k, 0) != sc.get(k, 0) for k in kw | ks)
    emit("l1_multiplicity_mismatch_groups", l1_mult_mismatch)
    l1_contract = (kw == ks and l1_mult_mismatch == 0 and len(gw) == len(gs))
    emit("l1_key_contract", "PASS" if l1_contract else "FAIL")

    l1_cols = l1_cells = 0; l1_abs = l1_rel = 0.0
    if l1_contract:
        CW, cw_co = canonicalize(gw, feats, "_k")
        CS, cs_co = canonicalize(gs, feats, "_k")
        emit("l1_cells_coerced_to_nan_whole", cw_co)
        emit("l1_cells_coerced_to_nan_split", cs_co)
        emit("l1_canonical_row_keys_equal", str(CW.index.equals(CS.index)))
        if CW.index.equals(CS.index):
            l1_cols, l1_cells, l1_abs, l1_rel = compare(CW, CS, feats, args.abs_tol, "l1", emit)
        else:
            l1_contract = False
            emit("l1_key_contract", "FAIL")
    emit("l1_features_compared", len(feats))
    emit("l1_features_differing", l1_cols)
    emit("l1_total_differing_cells", l1_cells)
    emit("l1_worst_max_abs_diff", f"{l1_abs:.6e}")
    emit("l1_worst_max_rel_diff", f"{l1_rel:.6e}")

    nw = W.loc[W["gene_id"].isna()].copy(); ns = S.loc[S["gene_id"].isna()].copy()
    emit("l2_rows_whole", len(nw)); emit("l2_rows_split", len(ns))
    vnw, vns = set(nw["_vh"]), set(ns["_vh"])
    emit("l2_variants_whole", len(vnw)); emit("l2_variants_split", len(vns))
    emit("l2_variants_only_in_whole", len(vnw - vns))
    emit("l2_variants_only_in_split", len(vns - vnw))
    wc2 = nw["_vh"].value_counts().to_dict(); sc2 = ns["_vh"].value_counts().to_dict()
    l2_mult_mismatch = sum(wc2.get(k, 0) != sc2.get(k, 0) for k in vnw | vns)
    emit("l2_multiplicity_mismatch_variants", l2_mult_mismatch)
    l2_contract = (vnw == vns and l2_mult_mismatch == 0 and len(nw) == len(ns))
    emit("l2_key_contract", "PASS" if l2_contract else "FAIL")

    l2_cols = l2_cells = 0; l2_abs = l2_rel = 0.0
    if l2_contract:
        NW, nw_co = canonicalize(nw, feats, "_vh")
        NS, ns_co = canonicalize(ns, feats, "_vh")
        emit("l2_cells_coerced_to_nan_whole", nw_co)
        emit("l2_cells_coerced_to_nan_split", ns_co)
        emit("l2_canonical_row_keys_equal", str(NW.index.equals(NS.index)))
        if NW.index.equals(NS.index):
            l2_cols, l2_cells, l2_abs, l2_rel = compare(NW, NS, feats, args.abs_tol, "l2", emit)
        else:
            l2_contract = False
            emit("l2_key_contract", "FAIL")
    emit("l2_features_compared", len(feats))
    emit("l2_features_differing", l2_cols)
    emit("l2_total_differing_cells", l2_cells)
    emit("l2_worst_max_abs_diff", f"{l2_abs:.6e}")
    emit("l2_worst_max_rel_diff", f"{l2_rel:.6e}")

    amb = 0
    if "Consequence_splice_region_variant" in nw.columns:
        for _, sub in nw.groupby("_vh", sort=False):
            if len(sub) > 1 and sub["Consequence_splice_region_variant"].nunique(dropna=False) > 1:
                amb += 1
    emit("within_arm_splice_region_ambiguous_variants", amb)
    emit("within_arm_note", "several candidate-transcript rows of one gene-unassigned variant already disagree on the splice-region flag inside a single arm; reported separately so it is never read as a sharding effect")

    if fatal or not l1_contract or not l2_contract:
        verdict = "FAIL_KEY_CONTRACT"
    elif l1_cols == 0 and l2_cols == 0:
        verdict = "PASS_SHARDING_EQUIVALENT"
    elif l1_cols == 0:
        verdict = "PASS_ANALYSIS_SPACE_ONLY_L2_DIFFERS"
    else:
        verdict = "FAIL_SHARDING_CHANGES_VALUES"

    emit("features_compared", len(feats))
    emit("features_differing", l1_cols + l2_cols)
    emit("total_differing_cells", l1_cells + l2_cells)
    emit("worst_max_abs_diff", f"{max(l1_abs, l2_abs):.6e}")
    emit("layer1_clean", str(l1_contract and l1_cols == 0))
    emit("layer2_clean", str(l2_contract and l2_cols == 0))
    emit("verdict", verdict)

    tmp = args.out + ".partial"
    with open(tmp, "w", encoding="utf-8") as fh:
        fh.write("key\tvalue\n" + "\n".join(lines) + "\n")
    os.replace(tmp, args.out)
    print(f"F2_V3 verdict={verdict} L1(contract={'PASS' if l1_contract else 'FAIL'} "
          f"diff_cols={l1_cols} cells={l1_cells}) L2(contract={'PASS' if l2_contract else 'FAIL'} "
          f"diff_cols={l2_cols} cells={l2_cells}) -> {args.out}")
    return 0 if verdict == "PASS_SHARDING_EQUIVALENT" else 1

if __name__ == "__main__":
    raise SystemExit(main())
