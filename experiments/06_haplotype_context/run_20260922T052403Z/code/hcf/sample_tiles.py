
import argparse
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

SEED = 20260922
MIN_MARKERS = 16
MAX_VARIANTS_PER_TILE = 4096

def _alloc(sizes, total):
    sizes = np.asarray(sizes, dtype=float)
    if sizes.sum() == 0:
        return np.zeros(len(sizes), dtype=int)
    raw = sizes / sizes.sum() * total
    base = np.floor(raw).astype(int)
    base = np.minimum(base, sizes.astype(int))
    rem = total - base.sum()
    order = np.argsort(-(raw - np.floor(raw)))
    i = 0
    while rem > 0 and i < len(order) * 4:
        j = order[i % len(order)]
        if base[j] < sizes[j]:
            base[j] += 1
            rem -= 1
        i += 1
    return base

def stratified_sample(df, total, seed, colname):
    rs = np.random.RandomState(seed)
    keys = sorted(df["stratum"].unique())
    sizes = [int((df["stratum"] == k).sum()) for k in keys]
    alloc = _alloc(sizes, total)
    chosen = []
    prob = pd.Series(0.0, index=df.index)
    for k, a, s in zip(keys, alloc, sizes):
        idx = df.index[df["stratum"] == k].values
        prob.loc[idx] = float(a) / s if s else 0.0
        if a > 0:
            chosen.extend(rs.choice(idx, size=int(a), replace=False).tolist())
    out = pd.Series(False, index=df.index)
    out.loc[chosen] = True
    return out, prob

def main(argv=None):
    ap = argparse.ArgumentParser("hcf.sample_tiles")
    ap.add_argument("--vardir", required=True)
    ap.add_argument("--out-parquet", required=True)
    ap.add_argument("--out-csv", required=True)
    ap.add_argument("--out-bench", required=True)
    ap.add_argument("--out-structure", required=True)
    ap.add_argument("--structure-max", type=int, default=300)
    ap.add_argument("--prediction-max", type=int, default=120)
    ap.add_argument("--bench-tiles", type=int, default=12)
    args = ap.parse_args(argv)

    frames = []
    for p in sorted(glob.glob(os.path.join(args.vardir, "chr*.tileinv.tsv"))):
        frames.append(pd.read_csv(p, sep="\t", dtype={"chrom": str}))
    df = pd.concat(frames, ignore_index=True)
    df["tile_end"] = df["tile_start"] + 100000 - 1
    df["tile_id"] = "chr" + df["chrom"].astype(str) + ":" + df["tile_start"].astype(str)
    df["build"] = "GRCh37"
    df["run_id"] = "HCF-20260922-v1/run_20260922T052403Z"
    df["protocol_version"] = "HCF-20260922-v1"

    reason = np.where(df["n_primary"] < MIN_MARKERS,
                      "TOO_FEW_PRIMARY_MARKERS(<%d)" % MIN_MARKERS,
                      np.where(df["n_primary"] > MAX_VARIANTS_PER_TILE,
                               "DENSE_WINDOW_RESOURCE_LIMIT(>%d)" % MAX_VARIANTS_PER_TILE, ""))
    df["exclusion_reason"] = reason
    df["eligible"] = df["exclusion_reason"] == ""

    el = df[df["eligible"]].copy()
    q = np.quantile(el["n_primary"].values, [0.25, 0.5, 0.75])
    el["variant_count_quartile"] = np.digitize(el["n_primary"].values, q, right=True)
    el["stratum"] = el["chrom"].astype(str) + "_q" + el["variant_count_quartile"].astype(str)

    sel_s, prob_s = stratified_sample(el, min(args.structure_max, len(el)), SEED, "structure")
    el["selected_structure"] = sel_s
    el["structure_inclusion_probability"] = prob_s
    st = el[el["selected_structure"]].copy()
    sel_p, prob_p = stratified_sample(st, min(args.prediction_max, len(st)), SEED, "prediction")
    el["selected_prediction"] = False
    el.loc[st.index[sel_p.values], "selected_prediction"] = True
    el["prediction_inclusion_probability"] = 0.0
    el.loc[st.index, "prediction_inclusion_probability"] = prob_p.values

    rs = np.random.RandomState(SEED)
    bench = []
    per = args.bench_tiles // 4
    for qi in range(4):
        pool = el.index[el["variant_count_quartile"] == qi].values
        rs.shuffle(pool)
        seen = set()
        pick = []
        for ix in pool:
            ch = el.loc[ix, "chrom"]
            if ch in seen:
                continue
            seen.add(ch)
            pick.append(ix)
            if len(pick) == per:
                break
        bench.extend(pick)
    el["selected_benchmark"] = False
    el.loc[bench, "selected_benchmark"] = True

    for c in ("variant_count_quartile", "stratum", "selected_structure", "selected_prediction",
              "selected_benchmark", "structure_inclusion_probability",
              "prediction_inclusion_probability"):
        df[c] = el[c] if c in ("stratum",) else None
    df = df.merge(el[["tile_id", "variant_count_quartile", "stratum", "selected_structure",
                      "selected_prediction", "selected_benchmark",
                      "structure_inclusion_probability", "prediction_inclusion_probability"]],
                  on="tile_id", how="left", suffixes=("_drop", ""))
    df = df[[c for c in df.columns if not c.endswith("_drop")]]
    df["sampling_seed"] = SEED
    df["sampling_strata"] = "chromosome x primary_variant_count_quartile"
    df["variant_filter_id"] = "noncoding_primary: not CDS, not splice_sensitive(+-8bp), biallelic SNV/short-indel<=50bp, polymorphic, (TYPED or imputed R2>=0.8), liftover-quarantine-free"

    df.to_parquet(args.out_parquet, index=False)
    df.to_csv(args.out_csv, index=False)
    el[el["selected_benchmark"]][["chrom", "tile_start", "n_primary",
                                  "variant_count_quartile"]].to_csv(args.out_bench, sep="\t",
                                                                    index=False)
    el[el["selected_structure"]].sort_values(["chrom", "tile_start"])[
        ["chrom", "tile_start", "n_primary", "variant_count_quartile", "stratum",
         "selected_prediction"]].to_csv(args.out_structure, sep="\t", index=False)

    summary = {
        "n_tiles_total": int(len(df)),
        "n_tiles_eligible": int(df["eligible"].sum()),
        "n_excluded_too_few": int((df["exclusion_reason"].str.startswith("TOO_FEW")).sum()),
        "n_excluded_dense": int((df["exclusion_reason"].str.startswith("DENSE")).sum()),
        "n_primary_variants_total": int(df["n_primary"].sum()),
        "n_keymap_records_total": int(df["n_keymap"].sum()),
        "n_coding_cds": int(df["n_coding_cds"].sum()),
        "n_splice_sensitive": int(df["n_splice_sensitive"].sum()),
        "n_quarantine": int(df["n_quarantine"].sum()),
        "n_lowq_imputed": int(df["n_lowq_imputed"].sum()),
        "n_monomorphic": int(df["n_monomorphic"].sum()),
        "n_not_biallelic_snv_indel": int(df["n_not_biallelic_snv_indel"].sum()),
        "n_primary_typed": int(df["n_primary_typed"].sum()),
        "quartile_breaks_n_primary": [float(x) for x in q],
        "n_structure_selected": int(df["selected_structure"].fillna(False).sum()),
        "n_prediction_selected": int(df["selected_prediction"].fillna(False).sum()),
        "n_benchmark_selected": int(df["selected_benchmark"].fillna(False).sum()),
        "seed": SEED,
    }
    print(json.dumps(summary, indent=2))
    return 0

if __name__ == "__main__":
    sys.exit(main())
