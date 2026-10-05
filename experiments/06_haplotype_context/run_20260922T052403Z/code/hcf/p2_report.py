
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ARM_ORDER = ["B_COV", "B_DS", "B_IND", "B_PAIR", "H_EXACT", "H_CLUSTER", "H_PHASE"]

def main(argv=None):
    ap = argparse.ArgumentParser("hcf.p2_report")
    ap.add_argument("--comp", required=True)
    ap.add_argument("--locked", required=True)
    ap.add_argument("--grid", required=True)
    ap.add_argument("--sim", default="")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)

    d = pd.read_csv(args.comp)
    lk = json.load(open(args.locked))
    ok = d[d["status"] == "OK"].copy()
    unstable = ok["unstable_cvmin_extrapolation"].astype(str).isin(["True", "true"]) \
        if "unstable_cvmin_extrapolation" in ok.columns else pd.Series(False, index=ok.index)
    out = {"n_unstable_cvmin_rows_excluded": int(unstable.sum()),
           "unstable_cvmin_rows": ok.loc[unstable, ["tile_id", "representation_id"]]
           .to_dict("records")}
    okc = ok[~unstable]

    rows = []
    for a in ARM_ORDER:
        s = ok[ok["representation_id"] == a]
        if not len(s):
            rows.append(dict(arm=a, n_tiles=0))
            continue
        g = s["gain_vs_BCOV_normalized"].astype(float)
        gc = okc[okc["representation_id"] == a]["gain_vs_BCOV_normalized_cvmin"].astype(float)
        rows.append(dict(
            arm=a, n_tiles=int(len(s)),
            r2_B_median=float(s["r2_B"].median()), r2_B_min=float(s["r2_B"].min()),
            r2_B_max=float(s["r2_B"].max()),
            gain_median=float(g.median()), gain_p10=float(g.quantile(0.1)),
            gain_p90=float(g.quantile(0.9)), gain_max=float(g.max()),
            n_gain_gt0=int((g > 0).sum()),
            gain_cvmin_median=float(gc.median()), n_gain_cvmin_gt0=int((gc > 0).sum()),
            n_shrunk_to_cov=int(s["shrunk_to_covariate_only"].astype(str).isin(
                ["True", "true", "1"]).sum()),
            median_df=float(s["df_effective"].median()),
            median_lambda=float(s["lambda_selected"].median()),
            median_alpha=float(s["alpha_selected"].median()),
            n_nonconv_candidates_median=float(s["n_nonconverged_candidates"].median())))
    arm_tab = pd.DataFrame(rows)
    arm_tab.to_csv(os.path.join(args.outdir, "p2_arm_summary.csv"), index=False)
    out["arm_summary"] = rows

    win = {}
    tied_all = 0
    for t, g in ok.groupby("tile_id"):
        gg = g[g["representation_id"] != "B_COV"]
        if not len(gg):
            continue
        mn = gg["sse_B"].min()
        w = gg.loc[gg["sse_B"] <= mn + 1e-9, "representation_id"].tolist()
        if len(w) == len(gg):
            tied_all += 1
        for a in w:
            win[a] = win.get(a, 0) + 1
    out["winner_counts"] = win
    out["n_tiles_all_arms_tied"] = tied_all

    base_arms = ["B_DS", "B_IND", "B_PAIR"]
    chal = ["H_EXACT", "H_CLUSTER", "H_PHASE"]
    piv = ok.pivot_table(index="tile_id", columns="representation_id", values="sse_B")
    have = [a for a in base_arms if a in piv.columns]
    bb = piv[have].min(axis=1)
    beats = {}
    for a in chal:
        if a not in piv.columns:
            beats[a] = None
            continue
        v = piv[a]
        m = v.notna() & bb.notna()
        beats[a] = dict(n_tiles_compared=int(m.sum()),
                        n_beats_best_baseline=int((v[m] < bb[m] - 1e-9).sum()),
                        n_ties=int((np.abs(v[m] - bb[m]) <= 1e-9).sum()))
    out["challenger_vs_best_baseline"] = beats
    pivc = okc.pivot_table(index="tile_id", columns="representation_id", values="sse_B_cvmin")
    havec = [a for a in base_arms if a in pivc.columns]
    bbc = pivc[havec].min(axis=1)
    beatsc = {}
    for a in chal:
        if a not in pivc.columns:
            beatsc[a] = None
            continue
        v = pivc[a]
        m = v.notna() & bbc.notna()
        beatsc[a] = dict(n_tiles_compared=int(m.sum()),
                         n_beats_best_baseline=int((v[m] < bbc[m] - 1e-9).sum()),
                         n_ties=int((np.abs(v[m] - bbc[m]) <= 1e-9).sum()))
    out["challenger_vs_best_baseline_cvmin"] = beatsc

    lt = lk["tiles"]
    out["locked"] = dict(
        n_tiles=len(lt),
        baseline_arm_counts=pd.Series([t["baseline_arm"] for t in lt]).value_counts().to_dict(),
        challenger_arm_counts=pd.Series([t["challenger_arm"] for t in lt]).value_counts().to_dict(),
        n_identical_predictions=int(sum(1 for t in lt if t.get("identical_predictions"))),
        n_ld_pruned=int(sum(1 for t in lt if t.get("ld_pruned"))),
        n_phase_not_identifiable=int(sum(
            1 for t in lt if t["arm_status"].get("H_PHASE") != "OK")))

    gr = pd.read_csv(args.grid)
    def lam_counts(df):
        return dict((str(a), dict((str(k), int(v)) for k, v in
                                  g["lam"].value_counts().items()))
                    for a, g in df.groupby("representation_id"))
    out["one_se_lambda_counts"] = lam_counts(
        gr[gr["is_one_se"].astype(str).isin(["True", "true"])])
    out["cv_min_lambda_counts"] = lam_counts(
        gr[gr["is_cv_min"].astype(str).isin(["True", "true"])])
    out["n_nonconverged_candidate_cells"] = int(
        (~gr["converged"].astype(str).isin(["True", "true"])).sum())

    with open(os.path.join(args.outdir, "p2_report_numbers.json"), "w") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1, default=str)

    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    dat, lbl = [], []
    for a in ARM_ORDER[1:]:
        s = ok[ok["representation_id"] == a]["gain_vs_BCOV_normalized"].astype(float).dropna()
        if len(s):
            dat.append(s.values)
            lbl.append(a.replace("_", "-"))
    if dat:
        bp = ax.boxplot(dat, labels=lbl, showfliers=True, widths=0.6)
        ax.axhline(0.0, color="0.3", lw=0.8, ls="--")
    ax.set_ylabel("normalized gain vs B-COV on B\n(SSE difference / (n_B · s²_A))")
    ax.set_title("HC-P2 · per-tile B-split gain by representation arm\n"
                 "[실측-신규] 120 prediction tiles · TCHL · family_aware=false")
    ax.tick_params(axis="x", rotation=20)
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "fig_p2_arm_gain_distribution.png"), dpi=170)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.6, 3.8))
    keys = [a for a in ARM_ORDER[1:] if a in win]
    ax.bar([k.replace("_", "-") for k in keys], [win[k] for k in keys], color="0.45")
    ax.set_ylabel("tiles where arm attains minimum SSE on B")
    ax.set_title("HC-P2 · per-tile winner counts (ties counted for every tied arm)\n"
                 "[실측-신규] tiles with all arms tied: %d" % tied_all)
    ax.tick_params(axis="x", rotation=20)
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "fig_p2_winner_counts.png"), dpi=170)
    plt.close(fig)

    if args.sim and os.path.exists(args.sim):
        sm = pd.read_csv(args.sim)
        ph = sm[(sm["scenario"] == "SYN-PHASE") & sm["delta"].notna()]
        if len(ph):
            fig, ax = plt.subplots(figsize=(6.6, 4.0))
            for t, g in ph.groupby("tile_id"):
                g = g.sort_values("delta")
                ax.plot(g["delta"].astype(float), g["rate"].astype(float), marker="o",
                        lw=1.0, label=t)
                ax.fill_between(g["delta"].astype(float), g["wilson_lo"].astype(float),
                                g["wilson_hi"].astype(float), alpha=0.12)
            ax.set_xscale("log")
            ax.axhline(0.8, color="0.3", ls="--", lw=0.8)
            ax.set_xlabel("δ (synthetic incremental signal fraction; not ΔR²)")
            ax.set_ylabel("rejection rate (one-sided α=0.025)")
            ax.set_ylim(-0.03, 1.03)
            ax.set_title("SYN-PHASE detection rate with Wilson CI\n[시뮬] 20 replicates per point")
            ax.legend(fontsize=6, ncol=2)
            fig.tight_layout()
            fig.savefig(os.path.join(args.outdir, "fig_p2_syn_phase_power.png"), dpi=170)
            plt.close(fig)
    print(json.dumps(out, ensure_ascii=False, indent=1, default=str)[:6000])
    return 0

if __name__ == "__main__":
    sys.exit(main())
