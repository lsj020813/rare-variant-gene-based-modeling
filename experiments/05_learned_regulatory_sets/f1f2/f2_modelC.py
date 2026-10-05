#!/usr/bin/env python
import os, sys, json, time, zlib
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1"); os.environ.setdefault("MKL_NUM_THREADS", "1"); os.environ["F_JOBS"] = "1"
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from f_common import *
import multiprocessing as mp

NAME = sys.argv[1]; SENS = len(sys.argv) > 2 and sys.argv[2] == "sens"
NSCREEN = 10; NREP = 20; NW = int(os.environ.get("F_WORKERS", "4"))
OD = f"{OUT}/{NAME}"; os.makedirs(OD, exist_ok=True)
t0 = time.time()
d = load(f"Cset_{NAME}")
X, cols, key, dom, chrom, forced = d["X"], d["cols"], d["key"], d["domain"], d["chr"], d["forced"]
n = len(X)
groups = {g: idx.values for g, idx in pd.Series(np.arange(n)).groupby(dom)}
doms = sorted(groups); dom_chr = {g: chrom[groups[g]][0] for g in doms}; dom_forced = {g: int(forced[groups[g]].max()) for g in doms}
log(NAME, "Model C domains", len(doms), "random", sum(1 for g in doms if dom_forced[g] == 0))

def dom_rng(g, salt=0):
    return np.random.default_rng(SEED + salt + zlib.crc32(g.encode()) % 1000000)

def dists(Xd, colmask=None):
    Xe = Xd if colmask is None else Xd[:, colmask]
    return eucl_dist_matrix(Xe), bio_dist_matrix(bio_parts(Xd, colmask))

def pass0(g):
    Xd = X[groups[g]]; De, Db = dists(Xd); k = min(MIN_SAMPLES, len(Xd) - 1)
    return np.partition(De, k, axis=1)[:, k], np.partition(Db, k, axis=1)[:, k]

if SENS:
    prim = json.load(open(f"{OUT}/ub_inub/f2C_chosen.json")); EPS = prim["eps_realized"]; CHOSEN = prim["chosen"]
else:
    with mp.get_context("fork").Pool(NW) as pool:
        P0 = pool.map(pass0, doms, chunksize=8)
    de10 = np.concatenate([a for a, b in P0]); db10 = np.concatenate([b for a, b in P0])
    EPS = dict(eucl={str(q): float(np.quantile(de10, q)) for q in EPS_Q_GRID}, bio={str(q): float(np.quantile(db10, q)) for q in EPS_Q_GRID})
    log("eps grid", EPS, round(time.time() - t0))

def fit(algo, sim, prm, Xd, De, Db, seed=0, colfrac=1.0):
    D = De if sim == "eucl" else Db
    if algo == "dbscan":
        eps = EPS[sim][str(prm)] * (np.sqrt(colfrac) if sim == "eucl" else 1.0)
        return run_cluster("dbscan", eps, D=D)
    if algo == "ward":
        return run_cluster("ward", prm, X=Xd)
    return run_cluster(algo, prm, D=D, seed=seed)

def pass1(g):
    idx = groups[g]; Xd = X[idx]; m = len(Xd); De, Db = dists(Xd); rr = dom_rng(g, 1)
    out = dict(domain=g, labels={}, screen=[])
    subs = [np.sort(rr.choice(m, int(0.8 * m), replace=False)) for _ in range(NSCREEN)]
    for algo, sim in COMBOS:
        grid = EPS_Q_GRID if algo == "dbscan" else K_GRID_C
        for prm in grid:
            lab = fit(algo, sim, prm, Xd, De, Db); ok, info = valid_partition(lab)
            out["labels"][(algo, sim, prm)] = lab.astype(np.int16)
            aris = [adjusted_rand_score(lab[s], fit(algo, sim, prm, Xd[s], De[np.ix_(s, s)], Db[np.ix_(s, s)])) for s in subs]
            out["screen"].append(dict(domain=g, algo=algo, similarity=sim, param=prm, valid=ok, ari_med=float(np.median(aris)), **info))
    return out

def select(S, dom_subset=None):
    if dom_subset is not None:
        S = S[S.domain.isin(dom_subset)]
    agg = S.groupby(["algo", "similarity", "param"]).agg(ari_med=("ari_med", "median"), frac_valid=("valid", "mean"), n_dom=("domain", "size"),
                                                          n_clusters_med=("n_clusters", "median"), largest_frac_med=("largest_frac", "median"), noise_frac_med=("noise_frac", "median")).reset_index()
    ch = {}
    for algo, sim in COMBOS:
        A = agg[(agg.algo == algo) & (agg.similarity == sim) & (agg.frac_valid >= 0.5)]
        ch[(algo, sim)] = float(A.sort_values(["ari_med", "param"], ascending=[False, True]).param.iloc[0]) if len(A) else None
    return ch, agg

if not SENS:
    with mp.get_context("fork").Pool(NW) as pool:
        P1 = []
        for i, r in enumerate(pool.imap_unordered(pass1, doms, chunksize=2)):
            P1.append(r)
            if (i + 1) % 100 == 0:
                log("pass1", i + 1, round(time.time() - t0))
    S = pd.DataFrame([row for r in P1 for row in r["screen"]])
    labels = {r["domain"]: r["labels"] for r in P1}
    CHOSEN, agg = select(S)
    agg["chosen"] = [CHOSEN[(a, s)] == p for a, s, p in zip(agg.algo, agg.similarity, agg.param)]
    agg.insert(0, "matrix", NAME); agg.insert(1, "model", "C"); agg.to_csv(f"{OD}/f2C_param_grid.csv", index=False)
    json.dump(dict(chosen={f"{a}|{s}": v for (a, s), v in CHOSEN.items()}, eps_realized=EPS), open(f"{OD}/f2C_chosen.json", "w"))
    log("chosen", CHOSEN, round(time.time() - t0))
    hold = []
    def norm(p, algo):
        return None if p is None else (float(p) if algo == "dbscan" else int(p))
    def holdout_rows(tag, rep, keep_doms):
        ch2, _ = select(S, keep_doms)
        for (algo, sim) in COMBOS:
            p0, p1 = CHOSEN[(algo, sim)], ch2[(algo, sim)]
            if p0 is None:
                continue
            same = (p1 is not None) and (norm(p1, algo) == norm(p0, algo))
            aris = []
            if not same and p1 is not None:
                for g in keep_doms:
                    a = labels[g][(algo, sim, norm(p0, algo))]; b = labels[g][(algo, sim, norm(p1, algo))]; aris.append(adjusted_rand_score(a, b))
            else:
                aris = [1.0] * len(keep_doms)
            hold.append(dict(matrix=NAME, model="C", algo=algo, similarity=sim, ptype=tag, rep=str(rep), param_full=p0, param_holdout=p1, selection_unchanged=same,
                             n_domains=len(keep_doms), ari_domain_median=float(np.median(aris)) if aris else np.nan, ari_domain_q1=float(np.quantile(aris, 0.25)) if aris else np.nan))
    for cc in sorted(set(dom_chr.values()), key=lambda s: int(s) if s.isdigit() else 99):
        holdout_rows("chr_holdout", cc, [g for g in doms if dom_chr[g] != cc])
    rr = np.random.default_rng(SEED + 5)
    for r in range(NREP):
        holdout_rows("dom_holdout20", r, list(rr.choice(doms, int(0.8 * len(doms)), replace=False)))
    pd.DataFrame(hold).to_csv(f"{OD}/f2C_selection_holdout.csv", index=False)
else:
    CHOSEN = {(k.split("|")[0], k.split("|")[1]): v for k, v in CHOSEN.items()}
    labels = None

def feature_mask(sim, rr):
    m = np.ones(88, bool)
    if sim == "eucl":
        m[rr.choice(88, 18, replace=False)] = False
    else:
        m[20 + rr.choice(60, 12, replace=False)] = False; m[80 + rr.choice(8, 2, replace=False)] = False
    return m

def pass2(g):
    idx = groups[g]; Xd = X[idx]; m = len(Xd); De, Db = dists(Xd); rr = dom_rng(g, 2)
    inv, cnt, first = dedup_rows(Xd); single = cnt == 1; rep = np.zeros(m, bool); rep[first] = True
    runs = []; cons = {}; sizes = []
    for algo, sim in COMBOS:
        prm = CHOSEN[(algo, sim)]
        if prm is None:
            continue
        prm = float(prm) if algo == "dbscan" else int(prm)
        ref = fit(algo, sim, prm, Xd, De, Db); okr, infor = valid_partition(ref)
        Pd = bio_parts(Xd)
        diag = dict(ari_ref_vs_ccre_class=float(adjusted_rand_score(ref, Pd["ccre"])), ari_ref_vs_rep_class=float(adjusted_rand_score(ref, Pd["rep"])),
                    ari_ref_vs_tf_absent=float(adjusted_rand_score(ref, (Pd["ntf"] == 0).astype(int))), ari_ref_vs_ccre_none=float(adjusted_rand_score(ref, (Pd["ccre"] == CCRE_NONE).astype(int))))
        plan = [("row80", r, np.sort(rr.choice(m, int(0.8 * m), replace=False)), None, 0) for r in range(NREP)]
        if not SENS:
            plan += [("feat80", r, np.arange(m), feature_mask(sim, rr), 0) for r in range(NREP)]
            if algo == "spectral":
                plan += [("seed", r, np.arange(m), None, r + 1) for r in range(NREP)]
        M = np.zeros((m, m), np.uint16); Pn = np.zeros((m, m), np.uint16)
        for ptype, r, s, cm, sd in plan:
            if cm is None:
                lab = fit(algo, sim, prm, Xd[s], De[np.ix_(s, s)], Db[np.ix_(s, s)], seed=sd)
            else:
                De2, Db2 = dists(Xd, cm); lab = fit(algo, sim, prm, Xd[:, cm], De2, Db2, seed=sd, colfrac=cm.mean())
            okp, infop = valid_partition(lab); a = ref[s]; met = agree(a, lab)
            ss, uu = single[s], rep[s]
            met["ari_excl_dup"] = float(adjusted_rand_score(a[ss], lab[ss])) if ss.sum() > 10 else np.nan
            met["ari_unique_rep"] = float(adjusted_rand_score(a[uu], lab[uu])) if uu.sum() > 10 else np.nan
            runs.append(dict(domain=g, algo=algo, similarity=sim, param=prm, ptype=ptype, rep=r, n=len(s), valid=okp, p_n_clusters=infop["n_clusters"], p_largest_frac=infop["largest_frac"], p_noise_frac=infop["noise_frac"],
                             ref_valid=okr, ref_n_clusters=infor["n_clusters"], ref_largest_frac=infor["largest_frac"], ref_noise_frac=infor["noise_frac"], **met))
            Pn[np.ix_(s, s)] += 1
            for c_ in np.unique(lab[lab >= 0]):
                ii = s[lab == c_]; M[np.ix_(ii, ii)] += 1
        comp = consensus_labels(M, Pn)
        cons[(algo, sim)] = comp.astype(np.int16)
        u, c = np.unique(comp, return_counts=True)
        sizes.append(dict(domain=g, algo=algo, similarity=sim, n=m, n_unique_vectors=int(len(first)), frac_singleton_rows=float(single.mean()), n_consensus_ge2=int((c >= 2).sum()), n_consensus_ge10=int((c >= 10).sum()),
                          largest_frac=float(c.max() / m), frac_rows_singleton_consensus=float((c == 1).sum() / m), frac_rows_in_ge10=float(c[c >= 10].sum() / m), size_median_ge2=float(np.median(c[c >= 2])) if (c >= 2).any() else np.nan,
                          ref_n_clusters=infor["n_clusters"], ref_largest_frac=infor["largest_frac"], ref_noise_frac=infor["noise_frac"], **diag))
    return dict(domain=g, runs=runs, cons=cons, sizes=sizes, keys=key[idx])

with mp.get_context("fork").Pool(NW) as pool:
    P2 = []
    for i, r in enumerate(pool.imap_unordered(pass2, doms, chunksize=2)):
        P2.append(r)
        if (i + 1) % 100 == 0:
            log("pass2", i + 1, round(time.time() - t0))
R = pd.DataFrame([row for r in P2 for row in r["runs"]]); R["forced"] = R.domain.map(dom_forced); R["chr"] = R.domain.map(dom_chr)
R.to_csv(f"{OD}/f2C_domain_runs.csv.gz", index=False)
SZ = pd.DataFrame([row for r in P2 for row in r["sizes"]]); SZ["forced"] = SZ.domain.map(dom_forced)
SZ.to_csv(f"{OD}/f2C_consensus_sizes_domains.csv", index=False)
if not SENS:
    ca = []
    for r in P2:
        for (algo, sim), comp in r["cons"].items():
            ca.append(pd.DataFrame(dict(key=r["keys"], domain=r["domain"], algo=algo, similarity=sim, consensus_cluster=comp)))
    pd.concat(ca, ignore_index=True).to_csv(f"{OD}/consensus_clusters.csv.gz", index=False)

stab = []
for algo, sim in COMBOS:
    prm = CHOSEN[(algo, sim)]
    if prm is None:
        stab.append(dict(matrix=NAME, model="C", algo=algo, similarity=sim, ptype="NONE_VALID")); continue
    G0 = R[(R.algo == algo) & (R.similarity == sim)]
    for scope, G in [("all1000", G0), ("rand701", G0[G0.forced == 0])]:
        pts = [(pt, G[G.ptype == pt]) for pt in G.ptype.unique()] + [("POOLED_partition_perturbations", G)]
        for pt, H in pts:
            s = summarize(H.ari); dm = H.groupby("domain").ari.median(); rm = H.groupby("rep").ari.median()
            stab.append(dict(matrix=NAME, model="C", scope=scope, algo=algo, similarity=sim, param=prm, ptype=pt, n_runs=len(H), n_domains=H.domain.nunique(),
                             ari_median=s["median"], ari_q1=s["q1"], ari_q3=s["q3"], ari_ci_lo=s["ci_lo"], ari_ci_hi=s["ci_hi"], grade=grade(s["median"]),
                             ari_domain_median_median=float(dm.median()), ari_domain_median_q1=float(dm.quantile(0.25)), ari_domain_median_q3=float(dm.quantile(0.75)), frac_domains_ari_ge07=float((dm >= 0.7).mean()), frac_domains_ari_lt04=float((dm < 0.4).mean()),
                             ari_runlevel_median=float(rm.median()), nmi_median=summarize(H.nmi)["median"], jaccard_median=summarize(H.jac)["median"],
                             ari_excl_dup_median=summarize(H.ari_excl_dup)["median"], ari_excl_dup_q1=summarize(H.ari_excl_dup)["q1"], ari_excl_dup_q3=summarize(H.ari_excl_dup)["q3"], grade_excl_dup=grade(summarize(H.ari_excl_dup)["median"]),
                             ari_unique_rep_median=summarize(H.ari_unique_rep)["median"], frac_runs_valid=float(H.valid.mean()), frac_domains_ref_valid=float(H.groupby("domain").ref_valid.first().mean()),
                             ref_n_clusters_median=float(H.groupby("domain").ref_n_clusters.first().median()), ref_largest_frac_median=float(H.groupby("domain").ref_largest_frac.first().median()),
                             ref_noise_frac_median=float(H.groupby("domain").ref_noise_frac.first().median()), eps_realized=EPS[sim][str(prm)] if algo == "dbscan" else np.nan))
stab = pd.DataFrame(stab)
P = stab[(stab.ptype == "POOLED_partition_perturbations") & (stab.scope == "all1000")].sort_values("ari_median", ascending=False)
best = P.iloc[0] if len(P) else None
J = dict(rule="F2 PASS <=> best Model C combo pooled ARI median >=0.7 AND ARI(excl duplicate rows) median >=0.4",
         best_combo=(f"{best.algo}|{best.similarity}|{best.param}" if best is not None else None), best_ari_median=(float(best.ari_median) if best is not None else None),
         best_ari_excl_dup_median=(float(best.ari_excl_dup_median) if best is not None else None),
         F2_PASS=bool(best is not None and best.ari_median >= 0.7 and best.ari_excl_dup_median >= 0.4),
         note="best-of-7 selection is optimistic; all combos reported")
stab["F2_PASS"] = J["F2_PASS"]; stab["best_combo"] = J["best_combo"]
stab.to_csv(f"{OD}/f2C_stability.csv", index=False)
cs = []
for (algo, sim), G in SZ.groupby(["algo", "similarity"]):
    for scope, H in [("all1000", G), ("rand701", G[G.forced == 0])]:
        cs.append(dict(matrix=NAME, model="C", algo=algo, similarity=sim, scope=scope, n_domains=len(H), n_rows=int(H.n.sum()), n_consensus_ge2_domain_median=float(H.n_consensus_ge2.median()),
                       n_consensus_ge10_domain_median=float(H.n_consensus_ge10.median()), largest_frac_domain_median=float(H.largest_frac.median()), frac_rows_singleton_consensus_median=float(H.frac_rows_singleton_consensus.median()),
                       frac_rows_in_ge10_median=float(H.frac_rows_in_ge10.median()), size_median_ge2_median=float(H.size_median_ge2.median()), frac_singleton_rows_median=float(H.frac_singleton_rows.median()),
                       n_consensus_ge2_total=int(H.n_consensus_ge2.sum()), n_consensus_ge10_total=int(H.n_consensus_ge10.sum()),
                       ari_ref_vs_ccre_class_median=float(H.ari_ref_vs_ccre_class.median()), ari_ref_vs_rep_class_median=float(H.ari_ref_vs_rep_class.median()),
                       ari_ref_vs_tf_absent_median=float(H.ari_ref_vs_tf_absent.median()), ari_ref_vs_ccre_none_median=float(H.ari_ref_vs_ccre_none.median())))
pd.DataFrame(cs).to_csv(f"{OD}/f2C_consensus_sizes.csv", index=False)
json.dump(J, open(f"{OD}/f2C_judgement.json", "w"), indent=1)
open(f"{OD}/f2C.done", "w").write(json.dumps(dict(seconds=time.time() - t0, judgement=J)))
log("F2 C DONE", NAME, J, round(time.time() - t0))
