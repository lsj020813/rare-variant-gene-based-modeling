
import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy import linalg as slinalg
from scipy import stats as sstats

from hcf import p2
from hcf import p2_run

SIM_ALPHAS = (0.0,)
SIM_LAMBDAS = p2.LAMBDAS
SIM_DELTAS = (1e-4, 1e-3, 1e-2)
ALPHA_ONE_SIDED = 0.025
PROJ_RIDGE = 1e-3

def wilson(k, n, z=1.959963984540054):
    if n == 0:
        return (float("nan"), float("nan"))
    ph = k / float(n)
    d = 1.0 + z * z / n
    c = (ph + z * z / (2 * n)) / d
    h = z * np.sqrt(ph * (1 - ph) / n + z * z / (4.0 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))

def risk_null_meta(ctx):
    CA = np.column_stack([np.ones(ctx.nA), ctx.cov[ctx.idxA].astype(np.float64)])
    gam = np.linalg.lstsq(CA, ctx.y[ctx.idxA].astype(np.float64), rcond=None)[0]
    resid = ctx.y[ctx.idxA].astype(np.float64) - CA @ gam
    return dict(sigma=float(resid.std(ddof=CA.shape[1])), gamma_cov=gam.tolist())

def risk_null(ctx, n_rep=1000, n_boot_rep=100, n_boot=2000, seed=p2.SEED):
    rng = np.random.RandomState(seed)
    nB = ctx.nB
    cov = ctx.cov[ctx.idxB].astype(np.float64)
    CA = np.column_stack([np.ones(ctx.nA), ctx.cov[ctx.idxA].astype(np.float64)])
    gam = np.linalg.lstsq(CA, ctx.y[ctx.idxA].astype(np.float64), rcond=None)[0]
    resid = ctx.y[ctx.idxA].astype(np.float64) - CA @ gam
    sigma = float(resid.std(ddof=CA.shape[1]))
    m = np.column_stack([np.ones(nB), cov]) @ gam
    a = ctx.cov[ctx.idxB, 4].astype(np.float64)
    a = 0.05 * sigma * (a - a.mean()) / a.std()
    ct = ctx.cov[ctx.idxB, 2].astype(np.float64)
    nc = ctx.cov[ctx.idxB, 3].astype(np.float64)
    scale = np.where(ct == 1, 1.0, np.where(nc == 1, 1.3, 0.7))
    rows = []
    for tag, sc in (("homoscedastic_normal", np.ones(nB)), ("cohort_heteroscedastic", scale)):
        rej = 0
        zs = []
        for _r in range(n_rep):
            eps = rng.normal(scale=sigma * sc, size=nB)
            d = (-4.0 * eps * a) / ctx.s2A
            mu, se = d.mean(), d.std(ddof=1) / np.sqrt(nB)
            z = mu / se
            zs.append(z)
            if sstats.norm.sf(z) < ALPHA_ONE_SIDED:
                rej += 1
        lo, hi = wilson(rej, n_rep)
        rows.append(dict(scenario="SYN-RISKNULL", test="paired_mean_z_one_sided",
                         residual_model=tag, tile_id="fixed_predictors", delta="",
                         n_replicates=n_rep, n_reject=rej, rate=rej / float(n_rep),
                         wilson_lo=lo, wilson_hi=hi,
                         mc_se=np.sqrt((rej / float(n_rep)) * (1 - rej / float(n_rep)) / n_rep),
                         nominal_alpha=ALPHA_ONE_SIDED, mean_z=float(np.mean(zs)),
                         sd_z=float(np.std(zs, ddof=1)), n_people_B=nB,
                         status="OK", limitation="fixed predictors; no model fitting"))
        rejb = 0
        for _r in range(n_boot_rep):
            eps = rng.normal(scale=sigma * sc, size=nB)
            d = (-4.0 * eps * a) / ctx.s2A
            W = rng.multinomial(nB, np.full(nB, 1.0 / nB), size=n_boot).astype(np.float32)
            bm = (W @ d) / nB
            if np.mean(bm <= 0.0) < ALPHA_ONE_SIDED:
                rejb += 1
        lo, hi = wilson(rejb, n_boot_rep)
        rows.append(dict(scenario="SYN-RISKNULL", test="paired_person_bootstrap_%d" % n_boot,
                         residual_model=tag, tile_id="fixed_predictors", delta="",
                         n_replicates=n_boot_rep, n_reject=rejb, rate=rejb / float(n_boot_rep),
                         wilson_lo=lo, wilson_hi=hi,
                         mc_se=np.sqrt((rejb / float(n_boot_rep)) *
                                       (1 - rejb / float(n_boot_rep)) / n_boot_rep),
                         nominal_alpha=ALPHA_ONE_SIDED, mean_z="", sd_z="", n_people_B=nB,
                         status="OK",
                         limitation="bootstrap replicates reduced to %d for budget" % n_boot_rep))
    return rows, dict(sigma=sigma, gamma_cov=gam.tolist())

def _xty_blocked(U, rows, gi, Y, block=8192):
    out = np.zeros((len(gi), Y.shape[1]), dtype=np.float64)
    for s in range(0, len(rows), block):
        idx = rows[s:s + block]
        Xb = np.asarray(U[np.ix_(idx, gi)], dtype=np.float64)
        out += Xb.T @ Y[idx]
    return out

def build_pass(ctx, blocks, train, tmpd, tag, arms):
    ZEX, _n1, _d = p2.state_block(blocks["keys_by_sw"], train)
    ZCL, _n2, _cst = p2_run.hapla_blocks(ctx, blocks["chrom"], blocks["tile_start"],
                                         blocks["pos"], train, tag, tmpd)
    if ZCL is None:
        ZCL = np.zeros((ctx.n, 0), dtype=np.float32)
    m_tile = blocks["DS"].shape[1]
    widths = [("NUIS", blocks["NUIS"].shape[1]), ("DS", m_tile), ("IND", 3 * m_tile),
              ("PAIR", blocks["PAIR"].shape[1]), ("ZEX", ZEX.shape[1]),
              ("ZCL", ZCL.shape[1]), ("Q", blocks["Q"].shape[1])]
    off, colmap = 0, {}
    for name, w in widths:
        colmap[name] = np.arange(off, off + w)
        off += w
    ncov = ctx.cov.shape[1]
    U = np.empty((ctx.n, 1 + ncov + off + 1), dtype=np.float32)
    U[:, 0] = 1.0
    U[:, 1:1 + ncov] = ctx.cov
    base = 1 + ncov
    for name, blk in (("NUIS", blocks["NUIS"]), ("DS", blocks["DS"]),
                      ("PAIR", blocks["PAIR"]), ("ZEX", ZEX), ("ZCL", ZCL),
                      ("Q", blocks["Q"])):
        if blk.shape[1]:
            U[:, base + colmap[name][0]:base + colmap[name][-1] + 1] = blk
    io = base + colmap["IND"][0]
    np.multiply(blocks["DS"], blocks["DS"], out=U[:, io:io + m_tile])
    U[:, io + m_tile:io + 2 * m_tile] = blocks["G"]
    np.multiply(blocks["G"], blocks["G"], out=U[:, io + 2 * m_tile:io + 3 * m_tile])
    U[:, -1] = 0.0
    ff = p2.FoldFit(U, ncov, train)
    sysm = {}
    for arm in arms:
        cols = (np.concatenate([colmap[k] for k in p2.ARM_BLOCKS[arm]])
                if p2.ARM_BLOCKS[arm] else np.zeros(0, dtype=np.int64))
        s = ff.arm_system(cols)
        sysm[arm] = None if s is None else dict(G=s[0], keep=s[2], sd=s[3])
    return dict(ff=ff, colmap=colmap, sysm=sysm)

def _solve_batch(ff, sysm, rows_fit, Y, rows_pred, lam):
    if sysm is None:
        C = np.asarray(ff.U[np.ix_(rows_fit, ff.c_idx)], dtype=np.float64)
        gam = ff.Acc_inv @ (C.T @ Y[rows_fit])
        return np.asarray(ff.U[np.ix_(rows_pred, ff.c_idx)], dtype=np.float64) @ gam
    gi = sysm["keep"] + 1 + ff.ncov
    D = 1.0 / sysm["sd"]
    Acx = ff.M[np.ix_(ff.c_idx, gi)]
    W = ff.Acc_inv @ Acx
    Cty = np.asarray(ff.U[np.ix_(rows_fit, ff.c_idx)], dtype=np.float64).T @ Y[rows_fit]
    Xty = _xty_blocked(ff.U, rows_fit, gi, Y)
    b = (Xty - W.T @ Cty) * D[:, None]
    A = np.array(sysm["G"], dtype=np.float64, order="C", copy=True)
    A.flat[::A.shape[0] + 1] += ff.n * lam
    cf = slinalg.cho_factor(A, lower=True, check_finite=False)
    del A
    Bs = slinalg.cho_solve(cf, b, check_finite=False)
    br = Bs * D[:, None]
    gam = ff.Acc_inv @ Cty - ff.Acc_inv @ (Acx @ br)
    Cp = np.asarray(ff.U[np.ix_(rows_pred, ff.c_idx)], dtype=np.float64)
    Xp = np.asarray(ff.U[np.ix_(rows_pred, gi)], dtype=np.float64)
    return Cp @ gam + Xp @ br

def cv_and_predict(ctx, blocks, arms, Ysyn, tmpd, tagpfx="sim", final_pass=None):
    R = Ysyn.shape[1]
    nl = len(SIM_LAMBDAS)
    fold_mse = dict((a, np.full((nl, p2.INNER_FOLDS, R), np.nan)) for a in arms)
    for fk in range(p2.INNER_FOLDS):
        val = ctx.folds[fk]
        train = np.setdiff1d(ctx.idxA, val)
        P = build_pass(ctx, blocks, train, tmpd, "%s_f%d" % (tagpfx, fk), arms)
        for arm in arms:
            for li, lam in enumerate(SIM_LAMBDAS):
                Yh = _solve_batch(P["ff"], P["sysm"][arm], train, Ysyn, val, lam)
                fold_mse[arm][li, fk, :] = np.mean((Ysyn[val] - Yh) ** 2, axis=0)
        del P
    sel = {}
    for arm in arms:
        cvm = np.nanmean(fold_mse[arm], axis=1)
        cvse = np.nanstd(fold_mse[arm], axis=1, ddof=1) / np.sqrt(p2.INNER_FOLDS)
        jmin = np.nanargmin(cvm, axis=0)
        thr = cvm[jmin, np.arange(R)] + cvse[jmin, np.arange(R)]
        pick = np.zeros(R, dtype=int)
        for r in range(R):
            ok = [i for i in range(nl) if cvm[i, r] <= thr[r]]
            pick[r] = max(ok, key=lambda i: SIM_LAMBDAS[i])
        sel[arm] = pick
    P = final_pass if final_pass is not None else build_pass(
        ctx, blocks, ctx.idxA, tmpd, "%s_fin" % tagpfx, arms)
    predB = dict((a, np.zeros((ctx.nB, R))) for a in arms)
    for arm in arms:
        for li, lam in enumerate(SIM_LAMBDAS):
            cols = np.where(sel[arm] == li)[0]
            if cols.size == 0:
                continue
            Yh = _solve_batch(P["ff"], P["sysm"][arm], ctx.idxA,
                              np.ascontiguousarray(Ysyn[:, cols]), ctx.idxB, lam)
            predB[arm][:, cols] = Yh
    lam_sel = dict((a, np.array([SIM_LAMBDAS[i] for i in sel[a]])) for a in arms)
    return predB, lam_sel

def residual_direction(ctx, P, v, arm):
    ff, sysm = P["ff"], P["sysm"][arm]
    v = np.asarray(v, dtype=np.float64).reshape(-1, 1)
    C = np.asarray(ff.U[:, ff.c_idx], dtype=np.float64)
    gam = ff.Acc_inv @ (C[ff.train].T @ v[ff.train])
    r = v - C @ gam
    if sysm is not None:
        gi = sysm["keep"] + 1 + ff.ncov
        D = 1.0 / sysm["sd"]
        Acx = ff.M[np.ix_(ff.c_idx, gi)]
        W = ff.Acc_inv @ Acx
        Xty = _xty_blocked(ff.U, ff.train, gi, r)
        b = (Xty - W.T @ (C[ff.train].T @ r[ff.train])) * D[:, None]
        A = np.array(sysm["G"], dtype=np.float64, order="C", copy=True)
        A.flat[::A.shape[0] + 1] += ff.n * PROJ_RIDGE
        beta = slinalg.cho_solve(slinalg.cho_factor(A, lower=True, check_finite=False), b,
                                 check_finite=False) * D[:, None]
        Xall = np.asarray(ff.U[:, gi], dtype=np.float64)
        r = r - ((Xall - C @ W) @ beta)
    r = r.ravel()
    sd = float(np.std(r[ctx.idxA]))
    return r, sd

def make_synthetic(ctx, r, delta, mu, sigma, n_rep, seed):
    rng = np.random.RandomState(seed)
    vr = float(np.var(r[ctx.idxA]))
    if vr <= 0 or not np.isfinite(vr):
        return None, 0.0
    base = float(np.var(mu[ctx.idxA])) + sigma ** 2
    gamma = np.sqrt(delta / (1.0 - delta) * base / vr)
    Y = np.empty((ctx.n, n_rep), dtype=np.float64)
    for j in range(n_rep):
        Y[:, j] = mu + gamma * r + rng.normal(scale=sigma, size=ctx.n)
    return Y, gamma

def _drop_marker(ctx, b, j):
    from hcf import tiles as htiles
    keep = [i for i in range(len(b["pos"])) if i != j]
    pos = b["pos"][keep]
    hap = np.ascontiguousarray(b["hap"][:, :, keep])
    DS = np.ascontiguousarray(b["DS"][:, keep])
    G = np.ascontiguousarray(b["G"][:, keep])
    var_keys = [b["var_keys"][i] for i in keep]
    remap = dict((o, i) for i, o in enumerate(keep))
    pairs = [(remap[a], remap[c], k) for (a, c, k) in b["pairs"]
             if a in remap and c in remap]
    PAIR = np.zeros((G.shape[0], len(pairs)), dtype=np.float32)
    Q = np.zeros((G.shape[0], len(pairs)), dtype=np.float32)
    for c, (a, d, _k) in enumerate(pairs):
        PAIR[:, c] = G[:, a] * G[:, d]
        Q[:, c] = (hap[:, 0, a].astype(np.float32) - hap[:, 1, a].astype(np.float32)) * \
                  (hap[:, 0, d].astype(np.float32) - hap[:, 1, d].astype(np.float32))
    subwins = htiles.subwindows(pos, markers=p2.MARKERS, step=p2.STEP, max_gap_bp=p2.MAX_GAP)
    keys_by_sw = p2.encode_subwindow_keys(hap, subwins, verify=False)
    out = dict(b)
    out.update(pos=pos, hap=hap, DS=DS, G=G, PAIR=PAIR, Q=Q, pairs=pairs,
               var_keys=var_keys, subwins=subwins, keys_by_sw=keys_by_sw)
    return out

def main(argv=None):
    ap = argparse.ArgumentParser("hcf.p2_sim")
    ap.add_argument("--run", required=True)
    ap.add_argument("--pheno", required=True)
    ap.add_argument("--tiles", required=True, help="대표 타일 TSV (chrom, tile_start)")
    ap.add_argument("--locked", required=True, help="locked_challengers.json")
    ap.add_argument("--out", required=True)
    ap.add_argument("--tmp", required=True)
    ap.add_argument("--hapla", default="")
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--risknull-reps", type=int, default=1000)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--skip-risknull", action="store_true")
    ap.add_argument("--no-tag", action="store_true")
    args = ap.parse_args(argv)
    p2_run.HAPLA = args.hapla or None
    args.tmp = os.path.join(args.tmp, "sim")
    os.makedirs(args.tmp, exist_ok=True)
    ctx = p2_run.Ctx(args)
    p2_run.log("sim ctx n=%d A=%d B=%d" % (ctx.n, ctx.nA, ctx.nB))

    if args.skip_risknull:
        meta = risk_null_meta(ctx)
    else:
        rows, meta = risk_null(ctx, n_rep=args.risknull_reps)
        pd.DataFrame(rows).to_csv(args.out, index=False)
        p2_run.log("risk-null done: %s" %
                   [(r["residual_model"], r["test"], r["rate"]) for r in rows])
    sigma = meta["sigma"]
    gam_cov = np.array(meta["gamma_cov"])
    mu_all = np.column_stack([np.ones(ctx.n), ctx.cov.astype(np.float64)]) @ gam_cov

    locked = json.load(open(args.locked))
    lk = dict((d["tile_id"], d) for d in locked["tiles"])
    tl = pd.read_csv(args.tiles, sep="\t")
    refits = 0
    for i in range(len(tl)):
        if (i % args.nshards) != args.shard:
            continue
        ch, ts = str(tl["chrom"].iloc[i]), int(tl["tile_start"].iloc[i])
        tile_id = "chr%s:%d" % (ch, ts)
        ent = lk.get(tile_id)
        if ent is None:
            continue
        base_arm, chal_arm = ent["baseline_arm"], ent["challenger_arm"]
        phase_arm = ("H_PHASE" if ent["arm_status"].get("H_PHASE") == "OK" else chal_arm)
        arms = sorted(set([base_arm, chal_arm, phase_arm]))
        tmpd = os.path.join(args.tmp, tile_id.replace(":", "_"))
        os.makedirs(tmpd, exist_ok=True)
        t0 = time.time()
        try:
            blocks, info = p2_run.build_tile_inputs(ctx, ch, ts, tile_id)
            if blocks is None:
                raise RuntimeError(info.get("status"))
            blocks["chrom"], blocks["tile_start"] = ch, ts
            scen_jobs = []
            from hcf import state as hstate
            cand = [(k, j) for j, k in enumerate(blocks["var_keys"])]
            afm = np.mean(blocks["G"][ctx.idxA], axis=0) / 2.0
            maf = np.minimum(afm, 1 - afm)
            cand = [(k, j) for (k, j) in cand if maf[j] >= 0.05]
            cand.sort(key=lambda t: hstate.key_hash(t[0]))
            causal_key, causal_j = (cand[0] if cand else (None, None))
            qidx = None
            if blocks["Q"].shape[1]:
                order = sorted(range(len(blocks["pairs"])),
                               key=lambda c: hstate.key_hash(blocks["pairs"][c][2]))
                qidx = order[0]
            Pfin = build_pass(ctx, blocks, ctx.idxA, tmpd, "fin", arms + ["B_COV"])
            if qidx is not None:
                r, sd = residual_direction(ctx, Pfin, blocks["Q"][:, qidx], base_arm)
                if sd <= 1e-10:
                    scen_jobs.append(("SYN-PHASE", None, None, "NO_RESIDUAL_DIRECTION"))
                else:
                    for d in SIM_DELTAS:
                        scen_jobs.append(("SYN-PHASE", d, r, "OK"))
            else:
                scen_jobs.append(("SYN-PHASE", None, None, "NO_SUPPORTED_PAIR"))
            if causal_j is not None:
                r, sd = residual_direction(ctx, Pfin, blocks["DS"][:, causal_j].astype(float),
                                           "B_COV")
                scen_jobs.append(("SYN-ADD", 1e-3, r, "OK" if sd > 1e-10
                                  else "NO_RESIDUAL_DIRECTION"))
                r_add = r
            else:
                scen_jobs.append(("SYN-ADD", None, None, "NO_COMMON_MARKER"))
            out_rows = []
            jobs_full = [j for j in scen_jobs if j[3] == "OK" and j[1] is not None]
            if jobs_full:
                Ys, gammas, tags = [], [], []
                for (sc, d, r, _st) in jobs_full:
                    Y, g = make_synthetic(ctx, r, d, mu_all, sigma, args.reps,
                                          p2.SEED + int(d * 1e6) + hash(sc) % 1000)
                    Ys.append(Y)
                    gammas.append(g)
                    tags.append((sc, d))
                Yall = np.concatenate(Ys, axis=1)
                predB, lamsel = cv_and_predict(ctx, blocks, arms, Yall, tmpd, "full",
                                               final_pass=Pfin)
                refits += 2 * Yall.shape[1]
                yB = ctx.y[ctx.idxB].astype(np.float64)
                for bi, (sc, d) in enumerate(tags):
                    sl = slice(bi * args.reps, (bi + 1) * args.reps)
                    ysyn_B = Yall[ctx.idxB, sl]
                    ca = phase_arm if sc == "SYN-PHASE" else chal_arm
                    m, se, z, pv = paired_test_syn(ysyn_B, predB[base_arm][:, sl],
                                                   predB[ca][:, sl], ctx.s2A)
                    out_rows.append(_row(sc, tile_id, ch, d, gammas[bi], m, pv, base_arm,
                                         ca, args.reps, ctx, lamsel, sl, info))
            del Pfin
            if causal_j is not None and not args.no_tag:
                r = r_add
                mblocks = _drop_marker(ctx, blocks, causal_j)
                mblocks["chrom"], mblocks["tile_start"] = ch, ts
                Y, g = make_synthetic(ctx, r, 1e-3, mu_all, sigma, args.reps, p2.SEED + 777)
                predB, lamsel = cv_and_predict(ctx, mblocks, arms, Y, tmpd, "tag")
                refits += 2 * Y.shape[1]
                m, se, z, pv = paired_test_syn(Y[ctx.idxB], predB[base_arm],
                                               predB[chal_arm], ctx.s2A)
                out_rows.append(_row("SYN-TAG", tile_id, ch, 1e-3, g, m, pv, base_arm,
                                     chal_arm, args.reps, ctx, lamsel,
                                     slice(0, args.reps), info,
                                     extra="causal marker masked from all arms"))
            for (sc, d, r, st) in scen_jobs:
                if st != "OK":
                    out_rows.append(dict(scenario=sc, tile_id=tile_id, chrom=ch, delta="",
                                         n_replicates=0, n_reject="", rate="", wilson_lo="",
                                         wilson_hi="", status=st,
                                         baseline_arm=base_arm, challenger_arm=chal_arm))
            pd.DataFrame(out_rows).to_csv(
                args.out.replace(".csv", ".%s.csv" % tile_id.replace(":", "_")), index=False)
            p2_run.log("%s sim done %.0fs refits=%d" % (tile_id, time.time() - t0, refits))
        except Exception as exc:
            import traceback
            p2_run.log("%s SIM FAILED %s: %s" % (tile_id, type(exc).__name__, exc))
            with open(os.path.join(args.tmp, "fail_%s.txt" % tile_id.replace(":", "_")),
                      "w") as fh:
                fh.write(traceback.format_exc())
    p2_run.log("total refits=%d (cap %d)" % (refits, 1200))
    return 0

def paired_test_syn(ysyn_B, yh_base, yh_chal, s2A):
    d = ((ysyn_B - yh_base) ** 2 - (ysyn_B - yh_chal) ** 2) / s2A
    m = d.mean(axis=0)
    se = d.std(axis=0, ddof=1) / np.sqrt(d.shape[0])
    z = m / np.maximum(se, 1e-300)
    return m, se, z, sstats.norm.sf(z)

def _row(sc, tile_id, ch, delta, gamma, m, pv, base_arm, chal_arm, reps, ctx, lamsel, sl,
         info, extra=""):
    k = int(np.sum(pv < ALPHA_ONE_SIDED))
    lo, hi = wilson(k, reps)
    st = "OK"
    if sc == "SYN-PHASE" and (hi < 0.8 or (lo < 0.8 < hi)):
        st = "POWER_UNRESOLVED" if lo < 0.8 < hi else "POWER_BELOW_TARGET"
    return dict(scenario=sc, tile_id=tile_id, chrom=ch, delta=delta, gamma=float(gamma),
                n_replicates=reps, n_reject=k, rate=k / float(reps), wilson_lo=lo,
                wilson_hi=hi, mc_se=np.sqrt((k / float(reps)) * (1 - k / float(reps)) / reps),
                mean_normalized_gain=float(np.mean(m)),
                sd_normalized_gain=float(np.std(m, ddof=1)),
                baseline_arm=base_arm, challenger_arm=chal_arm,
                lambda_base_median=float(np.median(lamsel[base_arm][sl])),
                lambda_chal_median=float(np.median(lamsel[chal_arm][sl])),
                n_people_B=ctx.nB, n_markers=info.get("n_markers"),
                n_supported_pairs=info.get("n_supported_pairs"),
                nominal_alpha=ALPHA_ONE_SIDED, tuning_grid="ridge_alpha0_7lambda_1SE",
                status=st, limitation=extra)

if __name__ == "__main__":
    sys.exit(main())

