#!/usr/bin/env python3

from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np

def make_design(n: int, rng: np.random.Generator) -> np.ndarray:
    age = rng.normal(55, 8, n)
    X = np.column_stack([
        np.ones(n),
        (age - age.mean()) / age.std(),
        ((age - age.mean()) / age.std()) ** 2,
        rng.integers(0, 2, n).astype(float),
        rng.normal(0, 1, (n, 10)),
    ])
    return X

def hc3_t(y: np.ndarray, X: np.ndarray, b: np.ndarray) -> float:
    D = np.column_stack([X, b])
    XtX_inv = np.linalg.pinv(D.T @ D)
    beta = XtX_inv @ (D.T @ y)
    resid = y - D @ beta
    h = np.einsum("ij,jk,ik->i", D, XtX_inv, D)
    h = np.clip(h, 0, 0.9999)
    omega = (resid / (1.0 - h)) ** 2
    meat = (D * omega[:, None]).T @ D
    V = XtX_inv @ meat @ XtX_inv
    se = np.sqrt(max(V[-1, -1], 1e-300))
    return abs(beta[-1] / se)

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=26070)
    ap.add_argument("--replicates", type=int, default=4000, help="null replicates per sparsity level")
    ap.add_argument("--bootstrap", type=int, default=2000, help="wild bootstrap draws per replicate")
    ap.add_argument("--modal-shares", default="0.985,0.9926,0.997,0.9995")
    ap.add_argument("--seed", type=int, default=20260726)
    ap.add_argument("--time-budget-s", type=int, default=5400)
    ap.add_argument("--out", default=os.environ.get("EXP_OUT"))
    args = ap.parse_args()
    if not args.out:
        sys.exit("no output path")

    rng = np.random.default_rng(args.seed)
    n = args.n
    X = make_design(n, rng)
    shares = [float(s) for s in args.modal_shares.split(",")]
    t0 = time.time()
    lines = []
    emit = lambda k, v: lines.append(f"{k}\t{v}")

    emit("n_samples", n)
    emit("covariates_plus_intercept", X.shape[1])
    emit("replicates_per_level", args.replicates)
    emit("bootstrap_draws", args.bootstrap)
    emit("seed", args.seed)

    from scipy import stats

    alphas = [1e-3, 1e-4, 2.63e-6]
    for share in shares:
        carriers = max(2, int(round(n * (1.0 - share))))
        wald_p = np.empty(args.replicates)
        boot_p = np.full(args.replicates, np.nan)
        n_boot_done = 0
        for r in range(args.replicates):
            b = np.zeros(n)
            idx = rng.choice(n, size=carriers, replace=False)
            b[idx] = rng.gamma(2.0, 1.0, carriers)
            y = X @ np.r_[0.3, np.zeros(X.shape[1] - 1)] + rng.normal(0, 1, n)
            t_obs = hc3_t(y, X, b)
            wald_p[r] = 2 * stats.norm.sf(t_obs)
            if time.time() - t0 < args.time_budget_s * 0.75 and r % 5 == 0:
                Xb = np.column_stack([X, b])
                XtX_inv = np.linalg.pinv(X.T @ X)
                bhat0 = XtX_inv @ (X.T @ y)
                e0 = y - X @ bhat0
                h0 = np.einsum("ij,jk,ik->i", X, XtX_inv, X)
                e0s = e0 / np.sqrt(np.clip(1.0 - h0, 1e-6, None))
                cnt = 0
                for _ in range(args.bootstrap):
                    ystar = X @ bhat0 + e0s * rng.choice([-1.0, 1.0], n)
                    if hc3_t(ystar, X, b) >= t_obs:
                        cnt += 1
                boot_p[r] = (1 + cnt) / (args.bootstrap + 1)
                n_boot_done += 1
            if time.time() - t0 > args.time_budget_s:
                wald_p = wald_p[: r + 1]
                boot_p = boot_p[: r + 1]
                emit(f"share_{share}_truncated_at_replicate", r + 1)
                break
        tag = f"share_{share}"
        emit(f"{tag}_carriers", carriers)
        emit(f"{tag}_replicates_done", int(wald_p.size))
        emit(f"{tag}_bootstrap_replicates_done", int(n_boot_done))
        for a in alphas:
            emit(f"{tag}_wald_tail_prob_at_{a:g}", f"{float((wald_p <= a).mean()):.6g}")
            emit(f"{tag}_wald_tail_ratio_vs_nominal_at_{a:g}",
                 f"{float((wald_p <= a).mean()) / a:.4f}")
        bp = boot_p[~np.isnan(boot_p)]
        if bp.size:
            for a in (1e-3, 1e-4):
                emit(f"{tag}_boot_tail_prob_at_{a:g}", f"{float((bp <= a).mean()):.6g}")
            emit(f"{tag}_boot_resolution_floor", f"{1.0/(args.bootstrap+1):.3g}")
        emit(f"{tag}_wald_p_min", f"{float(wald_p.min()):.3g}")
        emit(f"{tag}_elapsed_s", int(time.time() - t0))

    emit("verdict_note",
         "wald_tail_ratio_vs_nominal >> 1 means the asymptotic tail is anti-conservative "
         "at that alpha for that sparsity; ~1 means it is usable")
    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write("key\tvalue\n")
        fh.write("\n".join(lines) + "\n")
    print(f"G3_OK rows={len(lines)+1} elapsed_s={int(time.time()-t0)} out={args.out}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
