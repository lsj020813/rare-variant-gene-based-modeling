#!/usr/bin/env python
import os, sys, json, time, zlib
STAGE = sys.argv[2] if len(sys.argv) > 2 else "all"
_thr = "1" if STAGE == "C" else "4"
os.environ["OMP_NUM_THREADS"] = _thr; os.environ["OPENBLAS_NUM_THREADS"] = _thr; os.environ["MKL_NUM_THREADS"] = _thr
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from f_common import *
import multiprocessing as mp

NAME = sys.argv[1]
NREP = 20
PERDOM = int(os.environ.get("F_PERDOM", "100")); NQ = int(os.environ.get("F_NQ", "10000")); BLOCKNULL = os.environ.get("F_BLOCKNULL", "1") == "1"
OD = f"{OUT}/{NAME}"; os.makedirs(OD, exist_ok=True)
t0 = time.time()
d = load(f"Cset_{NAME}")
X, cols, key, dom, chrom, forced, mult = d["X"], d["cols"], d["key"], d["domain"], d["chr"], d["forced"], d["mult"]
n, p = X.shape
log(NAME, "loaded", X.shape)
dom_forced = pd.Series(forced).groupby(dom).max()
rand_dom = set(dom_forced.index[dom_forced == 0])
is_rand = np.array([x in rand_dom for x in dom])
log("domains", dom_forced.size, "random(non-forced)", len(rand_dom), "rows random", int(is_rand.sum()))

aud = []
def audit(Xs, keys, scope, vtype, label):
    inv, cnt, first = dedup_rows(Xs)
    m = len(Xs); nu = len(first)
    order = np.argsort(-cnt[first])
    g_sizes = cnt[first]
    df_ = pd.DataFrame(dict(g=inv, k=keys))
    gk = df_.groupby("g")["k"].nunique()
    rows_per_g = df_.groupby("g").size()
    same_key_dup = int((rows_per_g - gk).sum())
    cross_key_dup = int((gk - 1).sum())
    top = []
    for r in range(3):
        gi = inv[first[order[r]]] if r < nu else None
        if gi is None:
            break
        top.append(dict(size=int(g_sizes[order[r]]), frac=float(g_sizes[order[r]] / m), n_distinct_keys=int(gk.loc[gi]), desc=describe_vector(Xs[first[order[r]]], cols)))
    aud.append(dict(matrix=NAME, scope=scope, vector_type=vtype, label=label, n_rows=m, n_unique_vectors=nu, frac_unique=nu / m,
                    frac_rows_singleton=float((cnt == 1).mean()), n_groups_ge10=int((g_sizes >= 10).sum()), frac_rows_in_groups_ge10=float((cnt >= 10).mean()),
                    dup_rows_same_key=same_key_dup, dup_rows_cross_key=cross_key_dup,
                    largest_size=top[0]["size"], largest_frac=top[0]["frac"], largest_n_distinct_keys=top[0]["n_distinct_keys"], largest_desc=top[0]["desc"],
                    second_size=top[1]["size"] if len(top) > 1 else np.nan, second_desc=top[1]["desc"] if len(top) > 1 else "",
                    third_size=top[2]["size"] if len(top) > 2 else np.nan, third_desc=top[2]["desc"] if len(top) > 2 else ""))
    return inv, cnt

for scope, msk in [("all1000", np.ones(n, bool)), ("rand701", is_rand)]:
    audit(X[msk], key[msk], scope, "full88", "all 88 processed columns (rounded 1e-5)")
    audit(X[msk][:, :80], key[msk], scope, "discrete80", "cCRE/rep one-hot + cpg + 60 TF (annotation pattern only)")
    audit(X[msk][:, 20:80], key[msk], scope, "tf60", "60 TF indicators only")
pd.DataFrame(aud).to_csv(f"{OD}/duplicate_vector_audit.csv", index=False)
log("audit done", round(time.time() - t0))

def eff_rank(Xs):
    Xc = Xs.astype(np.float64) - Xs.mean(0)
    ev = np.linalg.eigvalsh(Xc.T @ Xc / (len(Xs) - 1))[::-1]
    ev = np.maximum(ev, 0); pr = ev / ev.sum()
    pz = pr[pr > 1e-12]
    H = -(pz * np.log(pz)).sum()
    cum = np.cumsum(pr)
    return dict(effective_rank_expH=float(np.exp(H)), participation_ratio=float(ev.sum() ** 2 / (ev ** 2).sum()), pc1_evr=float(pr[0]),
                n_comp_90=int(np.searchsorted(cum, 0.9) + 1), n_zero_eig=int((ev < 1e-10 * ev[0]).sum())), pr
er_rows = []; pca_curves = {}
for scope, msk in [("all1000", np.ones(n, bool)), ("rand701", is_rand)]:
    r, pr = eff_rank(X[msk]); r.update(matrix=NAME, scope=scope, n_rows=int(msk.sum())); er_rows.append(r); pca_curves[scope] = pr.tolist()
log("PCA done", er_rows[0])

rng = np.random.default_rng(SEED)
sub_idx = []
for g, idx in pd.Series(np.arange(n)).groupby(dom):
    idx = idx.values
    sub_idx.append(idx if len(idx) <= PERDOM else rng.choice(idx, PERDOM, replace=False))
sub_idx = np.sort(np.concatenate(sub_idx))
XU = X[sub_idx]; nU = len(XU); randU = is_rand[sub_idx]
log("Model U subsample", nU, "random rows", int(randU.sum()), "PERDOM", PERDOM, "NQ", NQ, "BLOCKNULL", BLOCKNULL)

def knn_stats(Xs, sim, k, q_idx=None):
    if sim == "eucl":
        Xr = Xs
    elif sim == "cos":
        Xr = Xs / (np.linalg.norm(Xs, axis=1, keepdims=True) + 1e-9)
    if sim in ("eucl", "cos"):
        if q_idx is None:
            dk, _ = knn_eucl(Xr, Xr, k, True)
        else:
            dk, _ = knn_eucl(Xr[q_idx], Xr, k + 1, False); dk = dk[:, 1:]
    else:
        P = bio_parts(Xs)
        if q_idx is None:
            dk, _ = knn_bio(P, P, k, True)
        else:
            dk, _ = knn_bio(sub_parts(P, q_idx), P, k + 1, False); dk = dk[:, 1:]
    return dk

def tend_metrics(dk, hop=None):
    d1 = dk[:, 0]; d10 = dk[:, min(9, dk.shape[1] - 1)]
    logdens = -np.log(d10 + 1e-9)
    bm = bimodality(logdens)
    out = dict(nn1_median=float(np.median(d1)), nn1_mean=float(d1.mean()), nn1_frac_zero=float((d1 <= 1e-6).mean()), nn10_median=float(np.median(d10)),
               dens_bc=bm["bc"], dens_kde_modes=bm["kde_modes"], hopkins=hop)
    return out, logdens

tend_rows = []; plot = dict(pca=pca_curves, nn1_hist={}, dens_hist={}, hopkins={})
def hist(v, bins):
    h, e = np.histogram(v, bins=bins); return dict(counts=h.tolist(), edges=e.tolist())
idrows = []

def run_modelU(Xs, tag, dedup_flag):
    m = len(Xs); rngU = np.random.default_rng(SEED + 7)
    lo, hi = Xs.min(0), Xs.max(0)
    for sim in ["eucl", "bio", "cos"]:
        k = 21
        dk = knn_stats(Xs, sim, k)
        hop = None
        if sim == "eucl":
            hops = [hopkins_eucl(Xs, 1000, rngU, lo, hi) for _ in range(NREP)]
            hop = float(np.median(hops)); plot["hopkins"][f"U_{tag}_real"] = hops
        real, ld_real = tend_metrics(dk, hop)
        idrows.append(dict(matrix=NAME, model="U", scope="all1000", similarity=sim, version=tag, kind="real", n=m,
                           twonn=twonn(dk[:, 0], dk[:, 1]), mle_k10=mle_id(dk, 10), mle_k20=mle_id(dk, 20), frac_nn1_zero=float((dk[:, 0] <= 1e-6).mean())))
        idrows.append(dict(matrix=NAME, model="U", scope="rand701", similarity=sim, version=tag, kind="real", n=int(randU_cur.sum()),
                           twonn=twonn(dk[randU_cur, 0], dk[randU_cur, 1]), mle_k10=mle_id(dk[randU_cur], 10), mle_k20=mle_id(dk[randU_cur], 20), frac_nn1_zero=float((dk[randU_cur, 0] <= 1e-6).mean())))
        nn1_bins = np.linspace(0, max(np.quantile(dk[:, 0], 0.995), 1e-3) * 1.5, 61)
        plot["nn1_hist"][f"U_{tag}_{sim}_real"] = hist(dk[:, 0], nn1_bins)
        dens_bins = np.linspace(np.quantile(ld_real, 0.002), np.quantile(ld_real, 0.998), 61)
        plot["dens_hist"][f"U_{tag}_{sim}_real"] = hist(ld_real, dens_bins)
        real701, _ = tend_metrics(dk[randU_cur], None)
        nulls = {"col": [], "block": []}
        for ntype in (["col", "block"] if (sim != "cos" and BLOCKNULL) else ["col"]):
            nrep = NREP
            for r in range(nrep):
                rr = np.random.default_rng(SEED + 1000 * (ntype == "block") + r)
                Y = shuffle_cols(Xs, rr, block=(ntype == "block"))
                q = rr.choice(m, min(NQ, m), replace=False)
                dkn = knn_stats(Y, sim, 21, q_idx=q)
                hopn = hopkins_eucl(Y, 1000, rr, lo, hi) if sim == "eucl" else None
                mets, ldn = tend_metrics(dkn, hopn)
                if r < 5:
                    idrows.append(dict(matrix=NAME, model="U", scope="all1000", similarity=sim, version=tag, kind=f"null_{ntype}_{r}", n=len(q),
                                       twonn=twonn(dkn[:, 0], dkn[:, 1]), mle_k10=mle_id(dkn, 10), mle_k20=mle_id(dkn, 20), frac_nn1_zero=float((dkn[:, 0] <= 1e-6).mean())))
                mets["ks_dens_vs_real"] = float(stats.ks_2samp(ld_real, ldn).statistic)
                nulls[ntype].append(mets)
                if r == 0:
                    plot["nn1_hist"][f"U_{tag}_{sim}_null{ntype}0"] = hist(dkn[:, 0], nn1_bins)
                    plot["dens_hist"][f"U_{tag}_{sim}_null{ntype}0"] = hist(ldn, dens_bins)
                if sim == "eucl":
                    plot["hopkins"].setdefault(f"U_{tag}_null{ntype}", []).append(hopn)
            log("U", tag, sim, ntype, "nulls done", round(time.time() - t0))
        for metric in ["hopkins", "nn1_median", "nn1_mean", "nn1_frac_zero", "nn10_median", "dens_bc", "dens_kde_modes"]:
            if real.get(metric) is None:
                continue
            for ntype, L in nulls.items():
                if not L:
                    continue
                nv = np.array([x[metric] for x in L], dtype=float)
                direction = "greater" if metric in ("hopkins", "dens_bc", "dens_kde_modes", "nn1_frac_zero") else "smaller"
                exceed = int((real[metric] > nv).sum()) if direction == "greater" else int((real[metric] < nv).sum())
                tend_rows.append(dict(matrix=NAME, model="U", scope="all1000", similarity=sim, version=tag, null_type=ntype, metric=metric, real=real[metric],
                                      real_rand701=real701.get(metric), null_median=float(np.median(nv)), null_min=float(nv.min()), null_max=float(nv.max()), n_null=len(nv),
                                      clustered_direction=direction, n_null_beaten=exceed, beats_all_nulls=bool(exceed == len(nv)),
                                      ks_dens_median=float(np.median([x["ks_dens_vs_real"] for x in L]))))
        log("U", tag, sim, "done", round(time.time() - t0), real)

if STAGE in ("U", "all"):
    randU_cur = randU
    run_modelU(XU, "withdup", False)
    invU, cntU, firstU = dedup_rows(XU)
    XUu = XU[np.sort(firstU)]; randU_cur = randU[np.sort(firstU)]
    log("Model U unique rows", len(XUu), "of", nU)
    run_modelU(XUu, "dedup", True)
    pd.DataFrame(tend_rows).to_csv(f"{OD}/cluster_tendency_modelU.csv", index=False)
    json.dump(dict(tend_rows=tend_rows, idrows=idrows, er_rows=er_rows, plot=plot, nU=int(nU), nU_unique=int(len(XUu))), open(f"{OD}/f1_stageU.json", "w"))
    log("stage U saved", round(time.time() - t0))
    if STAGE == "U":
        sys.exit(0)
else:
    st = json.load(open(f"{OD}/f1_stageU.json"))
    tend_rows, idrows, er_rows, plot = st["tend_rows"], st["idrows"], st["er_rows"], st["plot"]
    nU = st["nU"]; XUu = np.zeros((st["nU_unique"], 1))
    log("stage U loaded", len(tend_rows), len(idrows))

dom_groups = {g: idx.values for g, idx in pd.Series(np.arange(n)).groupby(dom)}

def domain_task(g):
    os.environ["OMP_NUM_THREADS"] = "1"; os.environ["OPENBLAS_NUM_THREADS"] = "1"
    idx = dom_groups[g]; Xd = X[idx]; m = len(idx)
    rr = np.random.default_rng(SEED + (zlib.crc32(g.encode()) % 1000000))
    res = dict(domain=g, n=m, chr=chrom[idx][0], forced=int(forced[idx].max()))
    inv, cnt, first = dedup_rows(Xd); Xdu = Xd[np.sort(first)]
    res["n_unique"] = len(Xdu); res["largest_dup_frac"] = float(cnt.max() / m)
    for tag, Xs in [("withdup", Xd), ("dedup", Xdu)]:
        ms = len(Xs)
        if ms < 12:
            continue
        lo, hi = Xs.min(0), Xs.max(0)
        De = eucl_dist_matrix(Xs); Db = bio_dist_matrix(bio_parts(Xs))
        def nnk(D, k):
            return np.sort(D + np.eye(ms, dtype=np.float32) * 1e9, axis=1)[:, :k]
        hmm = int(min(50, max(5, ms // 4)))
        hop_real = float(np.median([hopkins_eucl(Xs, hmm, rr, lo, hi) for _ in range(5)]))
        dke, dkb = nnk(De, 21), nnk(Db, 21)
        res[f"{tag}_hopkins"] = hop_real
        res[f"{tag}_eucl_nn1"] = float(np.median(dke[:, 0])); res[f"{tag}_bio_nn1"] = float(np.median(dkb[:, 0]))
        lde = -np.log(dke[:, min(9, ms - 2)] + 1e-9)
        if tag == "dedup":
            res["twonn_eucl"] = twonn(dke[:, 0], dke[:, 1]); res["twonn_bio"] = twonn(dkb[:, 0], dkb[:, 1])
            res["mle10_eucl"] = mle_id(dke, 10); res["mle10_bio"] = mle_id(dkb, 10)
            bm = bimodality(lde); res["dens_bc"] = bm["bc"]; res["dens_kde_modes"] = bm["kde_modes"]
        for ntype in ["col", "block"]:
            hn, e1, b1, ks = [], [], [], []
            for r in range(NREP):
                Y = shuffle_cols(Xs, rr, block=(ntype == "block"))
                hn.append(hopkins_eucl(Y, hmm, rr, lo, hi))
                Dn = eucl_dist_matrix(Y); dkn = nnk(Dn, 10)
                e1.append(float(np.median(dkn[:, 0])))
                ks.append(float(stats.ks_2samp(lde, -np.log(dkn[:, min(9, ms - 2)] + 1e-9)).statistic))
                if tag == "withdup" or ntype == "col":
                    Dbn = bio_dist_matrix(bio_parts(Y)); b1.append(float(np.median(np.sort(Dbn + np.eye(ms, dtype=np.float32) * 1e9, axis=1)[:, 0])))
            hn, e1 = np.array(hn), np.array(e1)
            res[f"{tag}_{ntype}_hopkins_null_med"] = float(np.median(hn)); res[f"{tag}_{ntype}_hopkins_nbeaten"] = int((hop_real > hn).sum())
            res[f"{tag}_{ntype}_eucl_nn1_null_med"] = float(np.median(e1)); res[f"{tag}_{ntype}_eucl_nn1_nbeaten"] = int((res[f"{tag}_eucl_nn1"] < e1).sum())
            res[f"{tag}_{ntype}_ks_dens_med"] = float(np.median(ks))
            if b1:
                b1 = np.array(b1); res[f"{tag}_{ntype}_bio_nn1_null_med"] = float(np.median(b1)); res[f"{tag}_{ntype}_bio_nn1_nbeaten"] = int((res[f"{tag}_bio_nn1"] < b1).sum())
    return res

log("Model C start", len(dom_groups))
with mp.get_context("fork").Pool(int(os.environ.get("F_WORKERS", "4"))) as pool:
    C = []
    for i, r in enumerate(pool.imap_unordered(domain_task, sorted(dom_groups.keys()), chunksize=4)):
        C.append(r)
        if (i + 1) % 100 == 0:
            log("Model C", i + 1, round(time.time() - t0))
C = pd.DataFrame(C).sort_values("domain"); C.to_csv(f"{OD}/f1_modelC_domains.csv", index=False)

for scope, msk in [("all1000", np.ones(len(C), bool)), ("rand701", (C.forced == 0).values)]:
    S = C[msk]
    for tag in ["withdup", "dedup"]:
        for ntype in ["col", "block"]:
            for metric, sim, direction in [("hopkins", "eucl", "greater"), ("eucl_nn1", "eucl", "smaller"), ("bio_nn1", "bio", "smaller")]:
                colr, coln, colb = f"{tag}_{metric}", f"{tag}_{ntype}_{metric}_null_med", f"{tag}_{ntype}_{metric}_nbeaten"
                if colb not in S:
                    continue
                ok = S[colb].notna()
                tend_rows.append(dict(matrix=NAME, model="C", scope=scope, similarity=sim, version=tag, null_type=ntype, metric=metric,
                                      real=float(S.loc[ok, colr].median()), real_rand701=np.nan, null_median=float(S.loc[ok, coln].median()), null_min=np.nan, null_max=np.nan, n_null=NREP,
                                      clustered_direction=direction, n_null_beaten=float(S.loc[ok, colb].median()), beats_all_nulls=float((S.loc[ok, colb] == NREP).mean()),
                                      ks_dens_median=float(S.loc[ok, f"{tag}_{ntype}_ks_dens_med"].median()), n_domains=int(ok.sum())))
    for sim in ["eucl", "bio"]:
        idrows.append(dict(matrix=NAME, model="C", scope=scope, similarity=sim, version="dedup", kind="real_domain_median", n=int(len(S)),
                           twonn=float(S[f"twonn_{sim}"].median()), mle_k10=float(S[f"mle10_{sim}"].median()), mle_k20=np.nan, frac_nn1_zero=np.nan,
                           twonn_q1=float(S[f"twonn_{sim}"].quantile(0.25)), twonn_q3=float(S[f"twonn_{sim}"].quantile(0.75))))
    tend_rows.append(dict(matrix=NAME, model="C", scope=scope, similarity="eucl", version="dedup", null_type="none", metric="dens_bc_frac_gt0.555", real=float((S.dens_bc > 0.555).mean()),
                          null_median=np.nan, n_null=0, clustered_direction="greater", n_domains=int(len(S))))
    tend_rows.append(dict(matrix=NAME, model="C", scope=scope, similarity="eucl", version="dedup", null_type="none", metric="dens_kde_modes_median", real=float(S.dens_kde_modes.median()),
                          null_median=np.nan, n_null=0, clustered_direction="greater", n_domains=int(len(S))))
    tend_rows.append(dict(matrix=NAME, model="C", scope=scope, similarity="eucl", version="withdup", null_type="none", metric="frac_unique_rows_domain_median", real=float((S.n_unique / S.n).median()),
                          null_median=np.nan, n_null=0, clustered_direction="none", n_domains=int(len(S))))
plot["modelC_hist"] = dict(hopkins_real=hist(C["withdup_hopkins"].dropna().values, np.linspace(0.5, 1, 51)), hopkins_null=hist(C["withdup_col_hopkins_null_med"].dropna().values, np.linspace(0.5, 1, 51)),
                           hopkins_real_dedup=hist(C["dedup_hopkins"].dropna().values, np.linspace(0.5, 1, 51)), hopkins_null_dedup=hist(C["dedup_col_hopkins_null_med"].dropna().values, np.linspace(0.5, 1, 51)),
                           twonn_eucl=hist(C["twonn_eucl"].dropna().values, np.linspace(0, 30, 61)), twonn_bio=hist(C["twonn_bio"].dropna().values, np.linspace(0, 30, 61)),
                           nn1_ratio=hist((C["dedup_eucl_nn1"] / C["dedup_col_eucl_nn1_null_med"]).dropna().values, np.linspace(0, 1.5, 61)))

T = pd.DataFrame(tend_rows)
def getv(model, version, metric, ntype="col", sim="eucl"):
    r = T[(T.model == model) & (T.version == version) & (T.metric == metric) & (T.null_type == ntype) & (T.similarity == sim) & (T.scope == "all1000")]
    return r.iloc[0] if len(r) else None
j = {}
for version in ["withdup", "dedup"]:
    h = getv("U", version, "hopkins"); e = getv("U", version, "nn1_median"); b = getv("U", version, "nn1_median", sim="bio")
    j[f"U_{version}_hopkins_beats_all"] = bool(h.beats_all_nulls) if h is not None else None
    j[f"U_{version}_nn1_eucl_beats_all"] = bool(e.beats_all_nulls) if e is not None else None
    j[f"U_{version}_nn1_bio_beats_all"] = bool(b.beats_all_nulls) if b is not None else None
j["F1_PASS_rule"] = "U withdup AND dedup: Hopkins beats all 20 col-nulls AND eucl 1-NN median below all 20 col-nulls"
j["F1_PASS"] = bool(all(j[k] for k in ["U_withdup_hopkins_beats_all", "U_withdup_nn1_eucl_beats_all", "U_dedup_hopkins_beats_all", "U_dedup_nn1_eucl_beats_all"]))
for version in ["withdup", "dedup"]:
    hc = getv("C", version, "hopkins"); ec = getv("C", version, "eucl_nn1")
    j[f"C_{version}_frac_domains_hopkins_beats_all"] = float(hc.beats_all_nulls) if hc is not None else None
    j[f"C_{version}_frac_domains_nn1_beats_all"] = float(ec.beats_all_nulls) if ec is not None else None
T["F1_PASS"] = j["F1_PASS"]
T.to_csv(f"{OD}/cluster_tendency.csv", index=False)
pd.DataFrame(idrows).to_csv(f"{OD}/intrinsic_dimension.csv", index=False)
pd.DataFrame(er_rows).to_csv(f"{OD}/effective_rank_recheck.csv", index=False)
plot["judgement"] = j; plot["n_rows"] = int(n); plot["nU"] = int(nU); plot["nU_unique"] = int(len(XUu)); plot["seconds"] = time.time() - t0
json.dump(plot, open(f"{OD}/f1_plotdata.json", "w"))
json.dump(j, open(f"{OD}/f1_judgement.json", "w"), indent=1)
open(f"{OD}/f1.done", "w").write(json.dumps(dict(seconds=time.time() - t0, judgement=j)))
log("F1 DONE", NAME, j)
