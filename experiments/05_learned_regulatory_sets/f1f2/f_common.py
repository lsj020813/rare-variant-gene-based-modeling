#!/usr/bin/env python
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import os, sys, json, gzip, time
import numpy as np
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import AgglomerativeClustering, DBSCAN, SpectralClustering
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from scipy import stats, sparse
from scipy.sparse.csgraph import connected_components

SEED = 20260922
PRIM = os.environ.get("F_PRIM", _config_path("${PROJECT_ROOT}/work/fset/primary"))
OUT = os.environ.get("F_OUT", _config_path("${PROJECT_ROOT}/work/fset/f1f2"))
W_CCRE, W_TF, W_REP, W_CPG, W_CONT = 2.0, 2.0, 0.5, 0.5, 1.0
CCRE_NONE, REP_NONE = 8, 9

def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)

def load(name):
    z = np.load(f"{PRIM}/{name}.npz", allow_pickle=True)
    cols = [str(c) for c in z["columns"]]
    check_cols(cols)
    d = dict(X=np.ascontiguousarray(z["X"], dtype=np.float32), cols=cols, key=z["key"].astype(str), domain=z["domain"].astype(str),
             chr=z["chr"].astype(str), forced=z["forced"].astype(int), mult=z["key_multiplicity"].astype(int))
    return d

def check_cols(cols):
    assert len(cols) == 88, len(cols)
    assert cols[0] == "ccre_PLS" and cols[8] == "ccre_none" and cols[9] == "rep_SINE" and cols[18] == "rep_none" and cols[19] == "cpg"
    assert all(c.startswith("tf_") for c in cols[20:80])
    assert cols[80:] == ["ccre_dist_rz", "rep_dist_rz", "cpg_dist_rz", "tf_n_rz", "map_k36_rz", "map_k36_missing", "gpn_score_rz", "gpn_score_missing"], cols[80:]

BLOCKS = [list(range(0, 9)) + [80], list(range(9, 19)) + [81], [19, 82]] + [[j] for j in range(20, 80)] + [[83], [84, 85], [86, 87]]

def shuffle_cols(X, rng, block=False):
    Y = X.copy()
    n = X.shape[0]
    if not block:
        for j in range(X.shape[1]):
            Y[:, j] = X[rng.permutation(n), j]
    else:
        for b in BLOCKS:
            p = rng.permutation(n)
            Y[:, b] = X[p][:, b]
    return Y

def bio_parts(X, colmask=None):
    P = dict(ccre=X[:, 0:9].argmax(1).astype(np.int16), rep=X[:, 9:19].argmax(1).astype(np.int16), cpg=(X[:, 19] > 0.5).astype(np.float32))
    tf = (X[:, 20:80] > 0.5).astype(np.float32)
    cont = np.clip(X[:, 80:88], -3, 3).astype(np.float32)
    if colmask is not None:
        tf = tf[:, colmask[20:80]]
        cont = cont[:, colmask[80:88]]
    P["tf"] = np.ascontiguousarray(tf); P["ntf"] = tf.sum(1)
    P["cont"] = np.ascontiguousarray(cont); P["cont2"] = (cont ** 2).sum(1); P["ncont"] = max(cont.shape[1], 1)
    return P

def sub_parts(P, idx):
    return dict(ccre=P["ccre"][idx], rep=P["rep"][idx], cpg=P["cpg"][idx], tf=P["tf"][idx], ntf=P["ntf"][idx], cont=P["cont"][idx], cont2=P["cont2"][idx], ncont=P["ncont"])

def bio_d2(A, B):
    inter = A["tf"] @ B["tf"].T
    union = A["ntf"][:, None] + B["ntf"][None, :]
    union -= inter
    np.maximum(union, 1.0, out=union)
    d2 = inter; d2 /= union
    d2 *= np.float32(-W_TF); d2 += np.float32(W_TF)
    c2 = A["cont"] @ B["cont"].T
    c2 *= np.float32(-2.0); c2 += A["cont2"][:, None]; c2 += B["cont2"][None, :]
    np.maximum(c2, 0, out=c2); c2 *= np.float32(W_CONT / (A["ncont"] * 4.0))
    d2 += c2; del c2
    same_c = (A["ccre"][:, None] == B["ccre"][None, :]); same_c &= (A["ccre"][:, None] != CCRE_NONE)
    same_c = ~same_c
    d2 += same_c; d2 += same_c
    same_r = (A["rep"][:, None] == B["rep"][None, :]); same_r &= (A["rep"][:, None] != REP_NONE)
    mism = (~same_r).view(np.int8)
    both = (A["cpg"][:, None] > 0.5) & (B["cpg"][None, :] > 0.5)
    mism += (~both).view(np.int8)
    d2 += np.multiply(mism, np.float32(W_REP), dtype=np.float32)
    return d2

def bio_dist_matrix(P):
    n = len(P["ccre"]); D = np.zeros((n, n), dtype=np.float32)
    step = 2000
    for i in range(0, n, step):
        D[i:i + step] = bio_d2(sub_parts(P, slice(i, i + step)), P)
    np.fill_diagonal(D, 0.0)
    return np.sqrt(np.maximum(D, 0))

def eucl_dist_matrix(X):
    X = X.astype(np.float32); s = (X ** 2).sum(1)
    D = s[:, None] + s[None, :] - 2.0 * (X @ X.T)
    np.fill_diagonal(D, 0.0)
    return np.sqrt(np.maximum(D, 0)).astype(np.float32)

def knn_eucl(Xq, Xr, k, self_included, step=2000):
    Xq = Xq.astype(np.float32); Xr = Xr.astype(np.float32)
    q2 = (Xq ** 2).sum(1); r2 = (Xr ** 2).sum(1); kk = k + (1 if self_included else 0)
    D = np.zeros((len(Xq), kk), np.float32); I = np.zeros((len(Xq), kk), np.int64)
    for s in range(0, len(Xq), step):
        d2 = Xq[s:s + step] @ Xr.T; d2 *= np.float32(-2.0); d2 += q2[s:s + step, None]; d2 += r2[None, :]
        if self_included:
            rows = np.arange(s, min(s + step, len(Xq))); d2[np.arange(len(rows)), rows] = -1.0
        idx = np.argpartition(d2, kk - 1, axis=1)[:, :kk]
        dd = np.take_along_axis(d2, idx, 1); o = np.argsort(dd, 1)
        D[s:s + step] = np.take_along_axis(dd, o, 1); I[s:s + step] = np.take_along_axis(idx, o, 1)
    if self_included:
        D, I = D[:, 1:], I[:, 1:]
    return np.sqrt(np.maximum(D, 0)), I

def knn_bio(Pq, Pr, k, self_included, step=1000):
    nq = len(Pq["ccre"]); kk = k + (1 if self_included else 0)
    D = np.zeros((nq, kk), dtype=np.float32); I = np.zeros((nq, kk), dtype=np.int64)
    for s in range(0, nq, step):
        d2 = bio_d2(sub_parts(Pq, slice(s, s + step)), Pr)
        if self_included:
            rows = np.arange(s, min(s + step, nq)); d2[np.arange(len(rows)), rows] = -1.0
        idx = np.argpartition(d2, kk - 1, axis=1)[:, :kk]
        dd = np.take_along_axis(d2, idx, 1); o = np.argsort(dd, 1)
        D[s:s + step] = np.take_along_axis(dd, o, 1); I[s:s + step] = np.take_along_axis(idx, o, 1)
    if self_included:
        D, I = D[:, 1:], I[:, 1:]
    return np.sqrt(np.maximum(D, 0)), I

def hopkins_eucl(X, m, rng, lo=None, hi=None):
    n = len(X)
    m = int(min(m, n // 2)) if n >= 4 else 1
    if lo is None:
        lo, hi = X.min(0), X.max(0)
    U = rng.uniform(lo, hi, size=(m, X.shape[1])).astype(np.float32)
    sidx = rng.choice(n, m, replace=False)
    nn = NearestNeighbors(n_neighbors=2, algorithm="brute", n_jobs=1).fit(X)
    u = nn.kneighbors(U, n_neighbors=1)[0][:, 0]
    w = nn.kneighbors(X[sidx], n_neighbors=2)[0][:, 1]
    return float(u.sum() / (u.sum() + w.sum() + 1e-12))

def twonn(d1, d2):
    ok = d1 > 0
    mu = d2[ok] / d1[ok]
    mu = np.sort(mu); N = len(mu)
    if N < 20:
        return np.nan
    keep = int(np.floor(N * 0.9))
    x = np.log(mu[:keep]); F = np.arange(1, keep + 1) / N
    y = -np.log(1 - F)
    return float((x * y).sum() / (x * x).sum()) if (x * x).sum() > 0 else np.nan

def mle_id(dk, k):
    T = dk[:, :k]
    ok = T[:, 0] > 0
    T = T[ok]
    if len(T) < 20:
        return np.nan
    lr = np.log(T[:, k - 1:k] / T[:, :k - 1])
    inv_m = lr.sum(1) / (k - 1)
    return float(1.0 / np.mean(inv_m)) if np.mean(inv_m) > 0 else np.nan

def bimodality(v):
    v = v[np.isfinite(v)]
    n = len(v)
    if n < 30:
        return dict(bc=np.nan, kde_modes=np.nan)
    g = stats.skew(v); kexc = stats.kurtosis(v)
    bc = (g ** 2 + 1) / (kexc + 3 * (n - 1) ** 2 / ((n - 2) * (n - 3)))
    z = (v - v.mean()) / (v.std() + 1e-12)
    try:
        kde = stats.gaussian_kde(z, bw_method="silverman")
        grid = np.linspace(z.min(), z.max(), 512); y = kde(grid)
        peaks = np.where((y[1:-1] > y[:-2]) & (y[1:-1] >= y[2:]) & (y[1:-1] > 0.05 * y.max()))[0]
        modes = int(len(peaks))
    except Exception:
        modes = np.nan
    return dict(bc=float(bc), kde_modes=modes)

def dedup_rows(X):
    Xr = np.ascontiguousarray(np.round(X, 5))
    v = Xr.view(np.dtype((np.void, Xr.dtype.itemsize * Xr.shape[1]))).ravel()
    _, first, inv, cnt = np.unique(v, return_index=True, return_inverse=True, return_counts=True)
    return inv.ravel(), cnt[inv.ravel()], first

def describe_vector(x, cols):
    parts = []
    L = len(x)
    if L == 60:
        tfs = [cols[20 + j][3:] for j in range(60) if x[j] > 0.5]
        return "TF{%d}=%s" % (len(tfs), ",".join(tfs) if tfs else "-")
    parts.append("cCRE=" + cols[int(np.argmax(x[0:9]))][5:])
    parts.append("rep=" + cols[9 + int(np.argmax(x[9:19]))][4:])
    parts.append("cpg=%d" % int(x[19] > 0.5))
    tfs = [cols[j][3:] for j in range(20, 80) if x[j] > 0.5]
    parts.append("TF{%d}=%s" % (len(tfs), ",".join(tfs) if tfs else "-"))
    if L >= 88:
        parts.append("cont(rz)=" + ",".join("%s:%.2f" % (cols[j].replace("_rz", ""), x[j]) for j in range(80, 88)))
    return " | ".join(parts)

K_GRID_U = [2, 3, 4, 5, 6, 8, 10, 12, 15, 20]
K_GRID_C = [2, 3, 4, 5, 6, 8, 10]
EPS_Q_GRID = [0.05, 0.10, 0.20, 0.35, 0.50]
MIN_SAMPLES = 10
COMBOS = [("ward", "eucl"), ("average", "eucl"), ("average", "bio"), ("dbscan", "eucl"), ("dbscan", "bio"), ("spectral", "eucl"), ("spectral", "bio")]

def knn_graph_affinity(D, k):
    n = D.shape[0]; k = min(k, n - 1)
    idx = np.argpartition(D, k, axis=1)[:, :k + 1]
    rows = np.repeat(np.arange(n), k + 1); cols = idx.ravel()
    m = rows != cols
    A = sparse.csr_matrix((np.ones(m.sum(), dtype=np.float32), (rows[m], cols[m])), shape=(n, n))
    A = A.maximum(A.T)
    return A

def run_cluster(algo, param, X=None, D=None, seed=0, n_init=3):
    n = X.shape[0] if X is not None else D.shape[0]
    if algo == "ward":
        k = int(min(param, n - 1))
        return AgglomerativeClustering(n_clusters=k, linkage="ward").fit(X).labels_
    if algo == "average":
        k = int(min(param, n - 1))
        import sklearn
        kw = {"metric": "precomputed"} if int(sklearn.__version__.split(".")[1]) >= 2 or int(sklearn.__version__.split(".")[0]) >= 2 else {"affinity": "precomputed"}
        return AgglomerativeClustering(n_clusters=k, linkage="average", **kw).fit(D.astype(np.float64)).labels_
    if algo == "dbscan":
        return DBSCAN(eps=float(param), min_samples=MIN_SAMPLES, metric="precomputed", n_jobs=1).fit(D.astype(np.float64)).labels_
    if algo == "spectral":
        k = int(min(param, n - 1))
        A = knn_graph_affinity(D, 15 if n > 1000 else 10)
        if n <= 1500:
            A = A.toarray()
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return SpectralClustering(n_clusters=k, affinity="precomputed", random_state=int(seed), assign_labels="kmeans", n_init=n_init).fit(A).labels_
    raise ValueError(algo)

def valid_partition(lab):
    m = lab >= 0
    if m.sum() == 0:
        return False, dict(n_clusters=0, largest_frac=1.0, noise_frac=1.0)
    u, c = np.unique(lab[m], return_counts=True)
    info = dict(n_clusters=int(len(u)), largest_frac=float(c.max() / len(lab)), noise_frac=float(1 - m.mean()))
    ok = len(u) >= 2 and info["largest_frac"] < 0.9 and info["noise_frac"] < 0.5
    return ok, info

def pair_jaccard(a, b):
    from sklearn.metrics.cluster import contingency_matrix
    C = contingency_matrix(a, b, sparse=True).astype(np.float64)
    nij = C.data
    A = (nij * (nij - 1) / 2).sum()
    ra = np.asarray(C.sum(1)).ravel(); cb = np.asarray(C.sum(0)).ravel()
    Ra = (ra * (ra - 1) / 2).sum(); Cb = (cb * (cb - 1) / 2).sum()
    den = Ra + Cb - A
    return float(A / den) if den > 0 else np.nan

def agree(a, b):
    return dict(ari=float(adjusted_rand_score(a, b)), nmi=float(normalized_mutual_info_score(a, b)), jac=pair_jaccard(a, b))

def eps_from_quantile(D, q):
    n = D.shape[0]; k = min(MIN_SAMPLES, n - 1)
    d10 = np.partition(D, k, axis=1)[:, k]
    return float(np.quantile(d10, q))

def summarize(vals):
    v = np.asarray([x for x in vals if x is not None and np.isfinite(x)], dtype=float)
    if len(v) == 0:
        return dict(n=0, median=np.nan, q1=np.nan, q3=np.nan, ci_lo=np.nan, ci_hi=np.nan, mean=np.nan)
    rng = np.random.default_rng(SEED)
    if len(v) >= 3:
        bs = np.median(rng.choice(v, size=(1000, len(v)), replace=True), axis=1)
        lo, hi = float(np.quantile(bs, 0.025)), float(np.quantile(bs, 0.975))
    else:
        lo = hi = np.nan
    return dict(n=int(len(v)), median=float(np.median(v)), q1=float(np.quantile(v, 0.25)), q3=float(np.quantile(v, 0.75)), ci_lo=lo, ci_hi=hi, mean=float(v.mean()))

def consensus_labels(M, Pn, thr=0.5):
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    co = M.astype(np.float32) / np.maximum(Pn, 1).astype(np.float32)
    co[Pn == 0] = 0.0
    D = 1.0 - co; np.fill_diagonal(D, 0.0); D = np.maximum(D, 0)
    if D.shape[0] < 2:
        return np.zeros(D.shape[0], np.int32)
    Z = linkage(squareform(D.astype(np.float64), checks=False), method="average")
    return (fcluster(Z, t=1.0 - thr + 1e-9, criterion="distance") - 1).astype(np.int32)

def grade(ari_median):
    if not np.isfinite(ari_median):
        return "NA"
    return "STRONG" if ari_median >= 0.7 else ("MODERATE" if ari_median >= 0.4 else "WEAK")
