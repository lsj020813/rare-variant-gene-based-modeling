
import argparse
import hashlib
import json
import os
import re
import sys

import numpy as np
import pandas as pd

SEED = 20260922
FRACTIONS = (0.6, 0.2, 0.2)

def make_split(sample_ids, seed=SEED):
    n = len(sample_ids)
    rs = np.random.RandomState(seed)
    perm = rs.permutation(n)
    nA = int(round(FRACTIONS[0] * n))
    nB = int(round(FRACTIONS[1] * n))
    lab = np.empty(n, dtype=object)
    lab[perm[:nA]] = "A"
    lab[perm[nA:nA + nB]] = "B"
    lab[perm[nA + nB:]] = "C"
    return lab

def batch_proxy(sample_ids):
    pref = [re.match(r"^[^0-9]*", s).group(0) for s in sample_ids]
    vc = pd.Series(pref).value_counts()
    usable = vc[vc >= 100]
    n = len(sample_ids)
    if 2 <= len(usable) <= 50:
        keep = set(usable.index)
        return np.array([p if p in keep else "OTHERPREFIX" for p in pref], dtype=object), \
            "sample_id_prefix (n_levels=%d)" % (len(usable) + 1)
    dec = (np.arange(n) * 10 // n)
    return np.array(["ORD%d" % d for d in dec], dtype=object), "vcf_sample_order_decile"

def main(argv=None):
    ap = argparse.ArgumentParser("hcf.split")
    ap.add_argument("--sample-list", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-parquet", default="")
    ap.add_argument("--out-meta", required=True)
    args = ap.parse_args(argv)

    ids = [l.strip() for l in open(args.sample_list) if l.strip()]
    lab = make_split(ids)
    bp, bp_kind = batch_proxy(ids)
    df = pd.DataFrame({"sample_index": np.arange(len(ids)), "sample_id": ids,
                       "split": lab, "batch_proxy": bp})
    df["family_group_id"] = ""
    df["family_aware"] = False
    df.to_csv(args.out_tsv, sep="\t", index=False)
    if args.out_parquet:
        df.to_parquet(args.out_parquet, index=False)
    meta = {
        "run_id": "HCF-20260922-v1/run_20260922T052403Z",
        "n_people": len(ids),
        "seed": SEED,
        "split_unit": "person",
        "family_aware": False,
        "fractions_requested": list(FRACTIONS),
        "counts": {k: int((lab == k).sum()) for k in ("A", "B", "C")},
        "sample_list_md5": hashlib.md5(open(args.sample_list, "rb").read()).hexdigest(),
        "batch_proxy_kind": bp_kind,
        "batch_proxy_levels": int(pd.Series(bp).nunique()),
        "blocked": ["BLOCKED_POPULATION_CONTROL"],
        "blocked_reason": ("no validated unrelated/kinship manifest; sparse GRM candidate "
                           "KINSHIP_CANDIDATE_FOUND_NOT_VALIDATED. Family/kinship-component "
                           "split not possible -> person-level split, EXPLORATORY_ONLY for any "
                           "final phenotype verdict."),
        "export_policy": "server-private; sample IDs and split membership are not exported",
        "used_by": ["HC-P2", "HC-V3"],
    }
    with open(args.out_meta, "w") as fh:
        json.dump(meta, fh, indent=2, ensure_ascii=False)
    print(json.dumps({k: v for k, v in meta.items() if k != "sample_list_md5"},
                     ensure_ascii=False))
    return 0

if __name__ == "__main__":
    sys.exit(main())
