#!/usr/bin/env python3

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

MAF_MAX = 1e-3
CADD_MIN = 5.0

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genes-ref", required=True)
    ap.add_argument("--annotations", required=True)
    ap.add_argument("--burden-genes-npy", required=False)
    ap.add_argument("--results", required=False)
    ap.add_argument("--out", default=os.environ.get("EXP_OUT"))
    args = ap.parse_args()
    if not args.out:
        sys.exit("no output path")

    lines = []
    emit = lambda k, v: lines.append(f"{k}\t{v}")

    gref = pd.read_parquet(args.genes_ref, columns=["id", "gene", "gene_type"])
    emit("reference_rows", len(gref))
    emit("reference_unique_gene_ids", int(gref["id"].nunique()))
    emit("reference_protein_coding",
         int((gref["gene_type"].astype(str) == "protein_coding").sum()))

    ann = pd.read_parquet(args.annotations, columns=["id", "gene_id", "MAF", "CADD_PHRED"])
    emit("annotation_rows", len(ann))
    emit("annotation_unique_variants", int(ann["id"].nunique()))
    genes_annotated = ann.loc[ann["gene_id"].notna(), "gene_id"].nunique()
    emit("genes_present_in_annotation", int(genes_annotated))

    maf = pd.to_numeric(ann["MAF"], errors="coerce")
    cadd = pd.to_numeric(ann["CADD_PHRED"], errors="coerce")
    elig = (maf > 0) & (maf < MAF_MAX) & (cadd > CADD_MIN) & ann["gene_id"].notna()
    genes_elig = ann.loc[elig, "gene_id"].nunique()
    emit("genes_with_eligible_variant", int(genes_elig))
    emit("eligible_variant_rows", int(elig.sum()))
    emit("eligible_unique_variants", int(ann.loc[elig, "id"].nunique()))

    if args.burden_genes_npy and os.path.exists(args.burden_genes_npy):
        gb = np.load(args.burden_genes_npy)
        emit("genes_on_burden_axis", int(np.unique(gb).size))
    if args.results and os.path.exists(args.results):
        res = pd.read_parquet(args.results)
        emit("genes_in_results", int(res["gene"].nunique()))

    for label, n in (("reference_all", int(gref["id"].nunique())),
                     ("genes_present_in_annotation", int(genes_annotated)),
                     ("genes_with_eligible_variant", int(genes_elig)),
                     ("declared_constant_19000", 19000)):
        emit(f"bonferroni_alpha_if_universe_{label}", f"{0.05/max(1,n):.6g}")

    emit("note", "chr21-only artifacts; genome-wide union requires the other chromosomes")

    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write("key\tvalue\n")
        fh.write("\n".join(lines) + "\n")
    print(f"F3B_OK rows={len(lines)+1} out={args.out}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
