
import argparse
import json
import sys

import numpy as np
import pandas as pd

from hcf import state as hstate

def main(argv=None):
    ap = argparse.ArgumentParser("hcf.p2_simtiles")
    ap.add_argument("--locked", required=True)
    ap.add_argument("--n", type=int, default=6)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    lk = json.load(open(args.locked))
    rows = []
    for t in lk["tiles"]:
        if t["arm_status"].get("H_PHASE") != "OK":
            continue
        rows.append(dict(tile_id=t["tile_id"], chrom=t["chrom"], tile_start=t["tile_start"],
                         n_supported_pairs=t["n_supported_pairs"], n_markers=t["n_markers"],
                         hsh=hstate.key_hash(t["tile_id"])))
    df = pd.DataFrame(rows)
    if not len(df):
        raise SystemExit("no tile with an identifiable phase arm")
    df["pair_stratum"] = pd.qcut(df["n_supported_pairs"].rank(method="first"), 3,
                                 labels=["p_low", "p_mid", "p_high"])
    med = df["n_markers"].median()
    df["size_stratum"] = np.where(df["n_markers"] >= med, "m_high", "m_low")
    df = df.sort_values("hsh")
    out, cells = [], sorted(set(zip(df["pair_stratum"].astype(str), df["size_stratum"])))
    while len(out) < args.n:
        added = False
        for cell in cells:
            if len(out) >= args.n:
                break
            sub = df[(df["pair_stratum"].astype(str) == cell[0]) &
                     (df["size_stratum"] == cell[1])]
            sub = sub[~sub["tile_id"].isin([o["tile_id"] for o in out])]
            if len(sub):
                out.append(sub.iloc[0].to_dict())
                added = True
        if not added:
            break
    o = pd.DataFrame(out)[["chrom", "tile_start", "tile_id", "n_supported_pairs", "n_markers",
                           "pair_stratum", "size_stratum"]]
    o.to_csv(args.out, sep="\t", index=False)
    print(o.to_string(index=False))
    return 0

if __name__ == "__main__":
    sys.exit(main())
