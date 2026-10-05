#!/usr/bin/env python
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1"); os.environ.setdefault("MKL_NUM_THREADS", "1"); os.environ["F_JOBS"] = "1"
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from f_common import *
import multiprocessing as mp

NAME = sys.argv[1]; SENS = len(sys.argv) > 2 and sys.argv[2] == "sens"
PER_DOM = 8; NSCREEN = 10; NREP = 20
OD = f"{OUT}/{NAME}"; os.makedirs(OD, exist_ok=True)
t0 = time.time()
d = load(f"Cset_{NAME}")
X, cols, key, dom, chrom, forced = d["X"], d["cols"], d["key"], d["domain"], d["chr"], d["forced"]
n = len(X)
dom_forced = pd.Series(forced).groupby(dom).max(); rand_dom = set(dom_forced.index[dom_forced == 0])
rng = np.random.default_rng(SEED)
sub = []
for g, idx in pd.Series(np.arange(n)).groupby(dom):
    idx = idx.values; sub.append(idx if len(idx) <= PER_DOM else rng.choice(idx, PER_DOM, replace=False))
sub = np.sort(np.concatenate(sub))
XU, keyU, domU, chrU = X[sub], key[sub], dom[sub], chrom[sub]
randU = np.array([g in rand_dom for g in domU]); nU = len(XU)
invU, cntU, firstU = dedup_rows(XU); singleton = cntU == 1; rep_mask = np.zeros(nU, bool); rep_mask[firstU] = True
log(NAME, "Model U subsample", nU, "unique", len(firstU), "singleton frac", round(singleton.mean(), 4))
PU = bio_parts(XU)
DE = eucl_dist_matrix(XU); DB = bio_dist_matrix(PU)
udoms = np.unique(domU); uchr = sorted(set(chrU), key=lambda s: int(s) if s.isdigit() else 99)

def dist_for(sim, idx, colmask=None):
    if colmask is None:
        return (DE if sim == "eucl" else DB)[np.ix_(idx, idx)]
    if sim == "eucl":
        return eucl_dist_matrix(XU[idx][:, colmask])
    return bio_dist_matrix(bio_parts(XU[idx], colmask))

def fit(algo, sim, param, idx, colmask=None, seed=0):
    D = dist_for(sim, idx, colmask)
    if algo == "dbscan":
        eps = eps_from_quantile(D, param); lab = run_cluster("dbscan", eps, D=D)
    elif algo == "ward":
        Xs = XU[idx] if colmask is None else XU[idx][:, colmask]; eps = np.nan; lab = run_cluster("ward", param, X=Xs)
    else:
        eps = np.nan; lab = run_cluster(algo, param, D=D, seed=seed)
    return lab, eps

def feature_mask(sim, rr):
    m = np.ones(88, bool)
    if sim == "eucl":
        m[rr.choice(88, 18, replace=False)] = False
    else:
        m[20 + rr.choice(60, 12, replace=False)] = False; m[80 + rr.choice(8, 2, replace=False)] = False
    return m

def compare(ref, lab, idx):
    r = ref[idx]; out = agree(r, lab)
    s = singleton[idx]; u = rep_mask[idx]; q = randU[idx]
    out["ari_excl_dup"] = float(adjusted_rand_score(r[s], lab[s])) if s.sum() > 10 else np.nan
    out["ari_unique_rep"] = float(adjusted_rand_score(r[u], lab[u])) if u.sum() > 10 else np.nan
    out["ari_rand701"] = float(adjusted_rand_score(r[q], lab[q])) if q.sum() > 10 else np.nan
    return out

def run_combo(combo):
    algo, sim = combo
    rr = np.random.default_rng(SEED + 17)
    grid = EPS_Q_GRID if algo == "dbscan" else K_GRID_U
    allidx = np.arange(nU); rows = []; screen = []
    if SENS:
        prim = json.load(open(f"{OUT}/ub_inub/f2U_chosen.json")); chosen = prim[f"{algo}|{sim}"]
        if chosen is None:
            return dict(combo=combo, screen=[], runs=[], chosen=None, consensus=None)
    else:
        refs = {}
        for prm in grid:
            lab, eps = fit(algo, sim, prm, allidx); ok, info = valid_partition(lab); refs[prm] = lab
            aris = []
            for r in range(NSCREEN):
                idx = np.sort(rr.choice(nU, int(0.8 * nU), replace=False))
                lb, _ = fit(algo, sim, prm, idx); aris.append(adjusted_rand_score(lab[idx], lb))
            screen.append(dict(matrix=NAME, model="U", algo=algo, similarity=sim, param=prm, eps_realized=eps, valid=ok, screen_ari_median=float(np.median(aris)),
                               screen_ari_q1=float(np.quantile(aris, 0.25)), screen_ari_q3=float(np.quantile(aris, 0.75)), **info))
        S = pd.DataFrame(screen); V = S[S.valid]
        chosen = float(V.sort_values(["screen_ari_median", "param"], ascending=[False, True]).param.iloc[0]) if len(V) else None
        if chosen is None:
            return dict(combo=combo, screen=screen, runs=[], chosen=None, consensus=None)
        if algo != "dbscan":
            chosen = int(chosen)
    ref, eps0 = fit(algo, sim, chosen, allidx); ok, info = valid_partition(ref)
    diag = dict(ari_ref_vs_ccre_class=float(adjusted_rand_score(ref, PU["ccre"])), ari_ref_vs_rep_class=float(adjusted_rand_score(ref, PU["rep"])),
                ari_ref_vs_tf_absent=float(adjusted_rand_score(ref, (PU["ntf"] == 0).astype(int))), ari_ref_vs_ccre_none=float(adjusted_rand_score(ref, (PU["ccre"] == CCRE_NONE).astype(int))))
    log("U", algo, sim, "chosen", chosen, info, diag, round(time.time() - t0))
    plan = []
    if SENS:
        for r in range(NREP):
            plan.append(("row80", r, np.sort(rr.choice(nU, int(0.8 * nU), replace=False)), None, 0))
    else:
        det = algo in ("ward", "average", "dbscan")
        for r in range(NREP if not det else 1):
            plan.append(("seed", r, allidx, None, r + 1))
        for r in range(NREP):
            plan.append(("row80", r, np.sort(rr.choice(nU, int(0.8 * nU), replace=False)), None, 0))
        for r in range(NREP):
            plan.append(("feat80", r, allidx, feature_mask(sim, rr), 0))
        for cc in uchr:
            plan.append((f"chr_holdout", cc, np.where(chrU != cc)[0], None, 0))
        for r in range(NREP):
            keep = rr.choice(udoms, int(0.8 * len(udoms)), replace=False); plan.append(("dom_holdout20", r, np.where(np.isin(domU, keep))[0], None, 0))
    M = np.zeros((nU, nU), np.uint16); Pn = np.zeros((nU, nU), np.uint16)
    for ptype, r, idx, cm, sd in plan:
        lab, eps = fit(algo, sim, chosen, idx, cm, sd)
        okp, infop = valid_partition(lab); met = compare(ref, lab, idx)
        deterministic = (ptype == "seed") and algo in ("ward", "average", "dbscan")
        rows.append(dict(matrix=NAME, model="U", algo=algo, similarity=sim, param=chosen, ptype=ptype, rep=str(r), n_rows=len(idx), eps_realized=eps, deterministic=deterministic,
                         valid=okp, **{("p_" + k): v for k, v in infop.items()}, **met))
        if not deterministic:
            Pn[np.ix_(idx, idx)] += 1
            for c_ in np.unique(lab[lab >= 0]):
                ii = idx[lab == c_]; M[np.ix_(ii, ii)] += 1
    comp = consensus_labels(M, Pn); ncomp = len(np.unique(comp))
    log("U", algo, sim, "battery done", len(plan), "consensus clusters", ncomp, round(time.time() - t0))
    return dict(combo=combo, screen=screen, runs=rows, chosen=chosen, ref=ref.astype(np.int32), consensus=comp.astype(np.int32), M=M, Pn=Pn, diag=diag)

with mp.get_context("fork").Pool(int(os.environ.get("F_WORKERS", "3"))) as pool:
    results = pool.map(run_combo, COMBOS, chunksize=1)

screen = pd.concat([pd.DataFrame(r["screen"]) for r in results if r["screen"]], ignore_index=True) if not SENS else pd.DataFrame()
runs = pd.concat([pd.DataFrame(r["runs"]) for r in results if r["runs"]], ignore_index=True)
if not SENS:
    screen.to_csv(f"{OD}/f2U_param_grid.csv", index=False)
    json.dump({f"{a}|{s}": r["chosen"] for r, (a, s) in zip(results, COMBOS)}, open(f"{OD}/f2U_chosen.json", "w"))
runs.to_csv(f"{OD}/f2U_runs.csv", index=False)

stab = []
for (algo, sim), r in zip(COMBOS, results):
    if r["chosen"] is None:
        stab.append(dict(matrix=NAME, model="U", algo=algo, similarity=sim, param=np.nan, ptype="NONE_VALID", n_runs=0)); continue
    R = runs[(runs.algo == algo) & (runs.similarity == sim)]
    refok, refinfo = valid_partition(r["ref"])
    groups = [(pt, R[R.ptype == pt]) for pt in R.ptype.unique()] + [("POOLED_nondeterministic", R[~R.deterministic])]
    for pt, G in groups:
        s = summarize(G.ari)
        stab.append(dict(matrix=NAME, model="U", algo=algo, similarity=sim, param=r["chosen"], ptype=pt, n_runs=len(G), deterministic=bool(G.deterministic.all()) if len(G) else False,
                         ari_median=s["median"], ari_q1=s["q1"], ari_q3=s["q3"], ari_ci_lo=s["ci_lo"], ari_ci_hi=s["ci_hi"], grade=grade(s["median"]),
                         nmi_median=summarize(G.nmi)["median"], jaccard_median=summarize(G.jac)["median"],
                         ari_excl_dup_median=summarize(G.ari_excl_dup)["median"], ari_excl_dup_q1=summarize(G.ari_excl_dup)["q1"], ari_excl_dup_q3=summarize(G.ari_excl_dup)["q3"],
                         ari_unique_rep_median=summarize(G.ari_unique_rep)["median"], ari_rand701_median=summarize(G.ari_rand701)["median"],
                         frac_runs_valid=float(G.valid.mean()) if len(G) else np.nan, ref_n_clusters=refinfo["n_clusters"], ref_largest_frac=refinfo["largest_frac"], ref_noise_frac=refinfo["noise_frac"],
                         eps_realized_ref=float(R.eps_realized.median()) if algo == "dbscan" else np.nan, n_rows_subsample=nU, frac_singleton_rows=float(singleton.mean()), **r["diag"]))
stab = pd.DataFrame(stab); stab.to_csv(f"{OD}/f2U_stability.csv", index=False)

cons_rows = []; assign = dict(key=keyU, domain=domU, singleton=singleton)
Mall = np.zeros((nU, nU), np.float64); Pall = np.zeros((nU, nU), np.float64)
for (algo, sim), r in zip(COMBOS, results):
    if r["chosen"] is None or SENS:
        continue
    comp = r["consensus"]; u, c = np.unique(comp, return_counts=True)
    assign[f"{algo}|{sim}"] = comp
    for scope, msk in [("all1000", np.ones(nU, bool)), ("rand701", randU), ("singleton_rows", singleton)]:
        uu, cc = np.unique(comp[msk], return_counts=True)
        cons_rows.append(dict(matrix=NAME, model="U", algo=algo, similarity=sim, scope=scope, n_rows=int(msk.sum()), n_consensus_clusters_ge2=int((cc >= 2).sum()), n_clusters_ge10=int((cc >= 10).sum()),
                              frac_rows_singleton_consensus=float((cc == 1).sum() / msk.sum()), largest_frac=float(cc.max() / msk.sum()), size_median_ge2=float(np.median(cc[cc >= 2])) if (cc >= 2).any() else np.nan,
                              size_q3_ge2=float(np.quantile(cc[cc >= 2], 0.75)) if (cc >= 2).any() else np.nan, frac_rows_in_ge10=float(cc[cc >= 10].sum() / msk.sum())))
    Mall += r["M"]; Pall += r["Pn"]
if not SENS:
    comp = consensus_labels(Mall, Pall); assign["POOLED_all_combos"] = comp.astype(np.int32)
    for scope, msk in [("all1000", np.ones(nU, bool)), ("rand701", randU), ("singleton_rows", singleton)]:
        uu, cc = np.unique(comp[msk], return_counts=True)
        cons_rows.append(dict(matrix=NAME, model="U", algo="POOLED", similarity="all", scope=scope, n_rows=int(msk.sum()), n_consensus_clusters_ge2=int((cc >= 2).sum()), n_clusters_ge10=int((cc >= 10).sum()),
                              frac_rows_singleton_consensus=float((cc == 1).sum() / msk.sum()), largest_frac=float(cc.max() / msk.sum()), size_median_ge2=float(np.median(cc[cc >= 2])) if (cc >= 2).any() else np.nan,
                              size_q3_ge2=float(np.quantile(cc[cc >= 2], 0.75)) if (cc >= 2).any() else np.nan, frac_rows_in_ge10=float(cc[cc >= 10].sum() / msk.sum())))
    np.savez_compressed(f"{OD}/f2U_consensus_assign.npz", **assign)
    pd.DataFrame(cons_rows).to_csv(f"{OD}/f2U_consensus_sizes.csv", index=False)
open(f"{OD}/f2U.done", "w").write(json.dumps(dict(seconds=time.time() - t0, nU=int(nU))))
log("F2 U DONE", NAME, round(time.time() - t0))
