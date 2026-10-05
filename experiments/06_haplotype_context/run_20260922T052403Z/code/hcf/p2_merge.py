
import argparse
import glob
import hashlib
import json
import os
import sys

import numpy as np
import pandas as pd

from hcf import p2

def pick(arms, sse, paired_se, present):
    av = [a for a in arms if a in present]
    if not av:
        return None, {}
    best = min(av, key=lambda a: sse[a])
    within = [a for a in av if (sse[a] - sse[best]) <= paired_se.get((a, best), 0.0)]
    if best not in within:
        within.append(best)
    chosen = min(within, key=lambda a: (p2.SIMPLICITY[a], sse[a]))
    return chosen, {"candidates": av, "sse": dict((a, float(sse[a])) for a in av),
                    "min_sse_arm": best, "within_1se": within}

def main(argv=None):
    ap = argparse.ArgumentParser("hcf.p2_merge")
    ap.add_argument("--run", required=True)
    ap.add_argument("--indir", required=True)
    ap.add_argument("--pheno", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args(argv)

    sp = pd.read_csv(os.path.join(args.run, "private/sample_split_manifest.private.tsv"),
                     sep="\t")
    ph = pd.read_csv(args.pheno, sep="\t")[["sample_id", "y"] + p2.COV_COLS]
    mg = sp[sp["split"].isin(["A", "B"])].merge(ph, on="sample_id", how="inner")
    mg = mg.dropna(subset=["y"] + p2.COV_COLS).sort_values("sample_index").reset_index(drop=True)
    yB = mg.loc[mg["split"] == "B", "y"].values.astype(np.float64)
    s2A = float(np.var(mg.loc[mg["split"] == "A", "y"].values.astype(np.float64), ddof=1))
    nB = yB.size

    bc = [pd.read_csv(f) for f in sorted(glob.glob(os.path.join(args.indir, "bc_*.csv")))]
    comp = pd.concat(bc, ignore_index=True) if bc else pd.DataFrame()
    grids, locks = [], []
    for f in sorted(glob.glob(os.path.join(args.indir, "lock_*.json"))):
        j = json.load(open(f))
        grids.extend(j.get("cv_grid", []))
        tile = j["tile_id"]
        npz = np.load(os.path.join(args.indir, "pred_%s.npz" % tile.replace(":", "_")))
        pred = dict((k, npz[k].astype(np.float64)) for k in npz.files)
        present = set(k for k in pred if not k.endswith("__cvmin"))
        sse = dict((a, float(np.sum((yB - pred[a]) ** 2))) for a in present)
        pse = {}
        for a in present:
            for b in present:
                if a == b:
                    continue
                d = (yB - pred[a]) ** 2 - (yB - pred[b]) ** 2
                pse[(a, b)] = float(d.std(ddof=1) / np.sqrt(nB)) * nB
        base_arm, base_dbg = pick(p2.BASELINE_ARMS, sse, pse, present)
        chal_arm, chal_dbg = pick(p2.CHALLENGER_ARMS, sse, pse, present)
        ent = dict(tile_id=tile, chrom=j["chrom"], tile_start=j["tile_start"],
                   baseline_arm=base_arm, challenger_arm=chal_arm,
                   baseline_selection=base_dbg, challenger_selection=chal_dbg,
                   arm_status=j["arm_status"], n_markers=j["info"].get("n_markers"),
                   n_supported_pairs=j["info"].get("n_supported_pairs"),
                   n_nuisance_cols=j["info"].get("n_nuisance_cols"),
                   ld_pruned=bool(j["info"].get("prune_applied")),
                   prune_r2=j["info"].get("prune_r2"),
                   n_design_cols_union=j["info"].get("n_design_cols_union"),
                   n_state_cols=j.get("n_zex_cols"), n_cluster_cols=j.get("n_zcl_cols"),
                   variant_keys_sha256=j["var_keys_sha256"],
                   nuisance_keys_sha256=j["nuis_keys_sha256"],
                   state_dictionary_sha256=j["state_dict_sha256"],
                   pair_keys=j["pair_keys"], hapla_status=j["hapla_note"])
        for role, arm in (("baseline", base_arm), ("challenger", chal_arm)):
            if arm is None:
                continue
            s = j["sel"][arm]
            fi = j["fitinfo"][arm]
            ent["%s_hyperparameters" % role] = dict(
                elastic_net_alpha=s["alpha"], lambda_=s["lam"],
                selection_rule="inner5foldCV_1SE_prefer_stronger_penalty",
                cv_mse_selected=s["cv_mse"], cv_mse_min=s["cv_mse_min"],
                cv_se_at_min=s["cv_se_at_min"],
                alpha_at_cv_min=s["alpha_min"], lambda_at_cv_min=s["lam_min"],
                n_nonconverged_candidates=s.get("n_nonconverged_candidates"))
            ent["%s_fit" % role] = dict(n_features=fi["n_features"], df_effective=fi["df"],
                                        n_nonzero=fi["nnz"],
                                        sse_B=sse[arm], r2_B=1.0 - sse[arm] /
                                        float(np.sum((yB - yB.mean()) ** 2)))
            ent["%s_blocks" % role] = list(p2.ARM_BLOCKS[arm])
        if base_arm and chal_arm:
            m, se, z, pv = p2.paired_gain(yB, pred[base_arm], pred[chal_arm], s2A)
            ent["B_paired_gain_challenger_minus_baseline"] = dict(
                normalized_mean=m, se=se, z=z, p_one_sided=pv,
                note="B 는 후보 고정용 탐색 표본이다. 이 값은 지지 판정이 아니다 (§8.4).")
            ent["identical_predictions"] = bool(np.allclose(pred[base_arm], pred[chal_arm]))
        locks.append(ent)

    os.makedirs(args.outdir, exist_ok=True)
    if len(comp):
        comp["unstable_cvmin_extrapolation"] = (
            comp["gain_vs_BCOV_normalized_cvmin"].abs() > 1.0)
    if len(comp):
        comp.to_csv(os.path.join(args.outdir, "baseline_comparison.csv"), index=False)
    if grids:
        pd.DataFrame(grids).to_csv(os.path.join(args.outdir, "baseline_cv_grid.csv"),
                                   index=False)
    payload = dict(
        run_id=p2.RUN_ID, protocol_version=p2.PROTOCOL_VERSION,
        evidence_class="[실측-신규]", stage="HC-P2 fit-development (A/B)",
        selection_scope="B only; C never opened in this track",
        split_manifest_sha256=p2.sha256_file(
            os.path.join(args.run, "private/sample_split_manifest.private.tsv")),
        config_sha256=p2.sha256_file(os.path.join(args.run, "config/resolved_config.yaml")),
        code_sha256=dict((os.path.basename(f), p2.sha256_file(f)) for f in sorted(
            glob.glob(os.path.join(args.run, "code/hcf/*.py")))),
        phenotype=dict(file="ref/pheno_v3/tchl_v3.tsv", column="y",
                       covariates=p2.COV_COLS, complete_case_n=int(len(mg)) ,
                       n_A=int((mg["split"] == "A").sum()), n_B=int(nB),
                       variance_in_A=s2A,
                       note="C split rows were never loaded by the fitting code"),
        family_aware=False,
        carried_blocked=["BLOCKED_POPULATION_CONTROL", "BLOCKED_PHASE_PROVENANCE_DOC_NOT_FOUND",
                         "BLOCKED_TYPED_QC_MANIFEST", "HDS_NOT_IN_MATCHED_SOURCE",
                         "KINSHIP_CANDIDATE_FOUND_NOT_VALIDATED", "MASTER_STATE_NOT_VERIFIED",
                         "ANNOTATION_RELEASE_NOT_RECORDED"],
        simplicity_rank=p2.SIMPLICITY,
        selection_rule=("per tile: min SSE on B among arms; any arm within one paired SE of "
                        "the minimum is treated as tied and the simpler arm (SIMPLICITY) is "
                        "locked. Rule fixed in code before A/B results were produced."),
        n_tiles=len(locks), tiles=locks)
    out = os.path.join(args.outdir, "locked_challengers.json")
    with open(out, "w") as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=1, default=str)
    h = p2.sha256_file(out)
    payload["self_sha256_of_previous_write"] = h
    with open(out, "w") as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=1, default=str)
    print("tiles=%d comparison_rows=%d grid_rows=%d sha=%s" % (len(locks), len(comp), len(grids), h))
    return 0

if __name__ == "__main__":
    sys.exit(main())
