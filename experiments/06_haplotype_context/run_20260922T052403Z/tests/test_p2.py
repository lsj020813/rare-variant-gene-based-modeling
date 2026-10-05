
import sys

import numpy as np

from hcf import p2
from hcf import phase as hphase
from hcf import state as hstate

FAILS = []

def check(name, cond, extra=""):
    if cond:
        print("PASS %s %s" % (name, extra))
    else:
        print("FAIL %s %s" % (name, extra))
        FAILS.append(name)

def _toy(n=400, p=12, seed=7):
    rng = np.random.RandomState(seed)
    X = rng.normal(size=(n, p)).astype(np.float32)
    C = np.column_stack([rng.normal(size=n), rng.randint(0, 2, n)]).astype(np.float32)
    y = (0.7 * X[:, 0] - 0.4 * X[:, 3] + 0.9 * C[:, 0] + rng.normal(scale=0.5, size=n)
         ).astype(np.float32)
    U = np.empty((n, 1 + C.shape[1] + p + 1), dtype=np.float32)
    U[:, 0] = 1.0
    U[:, 1:1 + C.shape[1]] = C
    U[:, 1 + C.shape[1]:-1] = X
    U[:, -1] = y
    return U, C, X, y

def t_ridge_matches_direct():
    U, C, X, y = _toy()
    n = X.shape[0]
    ff = p2.FoldFit(U, C.shape[1], np.arange(n))
    lam = 0.05
    fit = ff.solve(np.arange(X.shape[1]), 0.0, lam, want_df=True)
    yh = ff.predict(fit, np.arange(n))
    Cd = np.column_stack([np.ones(n), C]).astype(np.float64)
    sd = X.astype(np.float64).std(axis=0)
    Xs = X.astype(np.float64) / sd
    Z = np.column_stack([Cd, Xs])
    P = np.zeros((Z.shape[1], Z.shape[1]))
    P[Cd.shape[1]:, Cd.shape[1]:] = np.eye(Xs.shape[1]) * n * lam
    beta = np.linalg.solve(Z.T @ Z + P, Z.T @ y.astype(np.float64))
    yh2 = Z @ beta
    check("ridge_vs_direct_unpenalised_covariates",
          np.allclose(yh, yh2, atol=1e-6), "max|d|=%.2e" % np.max(np.abs(yh - yh2)))

def t_enet_matches_sklearn():
    from sklearn.linear_model import ElasticNet
    U, C, X, y = _toy(seed=11)
    n, p = X.shape
    sd = X.astype(np.float64).std(axis=0)
    mu = X.astype(np.float64).mean(axis=0)
    Xs = (X.astype(np.float64) - mu) / sd
    yc = y.astype(np.float64) - y.astype(np.float64).mean()
    G = Xs.T @ Xs
    b = Xs.T @ yc
    for alpha, lam in ((1.0, 0.02), (0.5, 0.05)):
        beta, _ni, _cv = p2.enet_gram_solve(G, b, yc, n, alpha, lam)
        en = ElasticNet(alpha=lam, l1_ratio=alpha, fit_intercept=False, max_iter=20000,
                        tol=1e-10)
        en.fit(Xs, yc)
        check("enet_gram_vs_sklearn_a%.1f_l%.2f" % (alpha, lam),
              np.allclose(beta, en.coef_, atol=1e-4),
              "max|d|=%.2e nnz=%d/%d" % (np.max(np.abs(beta - en.coef_)),
                                         int((beta != 0).sum()), int((en.coef_ != 0).sum())))

def t_cis_trans():
    hap = np.array([[[1, 1], [0, 0]], [[1, 0], [0, 1]]], dtype=np.int8)
    g0 = hap[0, 0] + hap[0, 1]
    g1 = hap[1, 0] + hap[1, 1]
    q = hphase.phase_contrast_q(hap, 0, 1)
    check("cis_trans_same_G_opposite_Q",
          np.array_equal(g0, g1) and q[0] == 1 and q[1] == -1,
          "G=%s Q=%s" % (g0.tolist(), q.tolist()))

def t_copy_swap_invariance():
    rng = np.random.RandomState(3)
    n, m = 300, 16
    hap = rng.randint(0, 2, size=(n, 2, m)).astype(np.int8)
    swapped = hap[:, ::-1, :].copy()
    k1 = p2.encode_subwindow_keys(hap, [(0, m)], verify=False)
    k2 = p2.encode_subwindow_keys(swapped, [(0, m)], verify=False)
    Z1, n1, _ = p2.state_block(k1, np.arange(n))
    Z2, n2, _ = p2.state_block(k2, np.arange(n))
    G1 = hap[:, 0, :] + hap[:, 1, :]
    G2 = swapped[:, 0, :] + swapped[:, 1, :]
    q1 = (hap[:, 0, 0] - hap[:, 1, 0]) * (hap[:, 0, 5] - hap[:, 1, 5])
    q2 = (swapped[:, 0, 0] - swapped[:, 1, 0]) * (swapped[:, 0, 5] - swapped[:, 1, 5])
    check("copy_swap_invariance_Z_G_Q",
          n1 == n2 and np.array_equal(Z1, Z2) and np.array_equal(G1, G2)
          and np.array_equal(q1, q2))

def t_single_marker_window():
    rng = np.random.RandomState(5)
    n = 200
    hap = rng.randint(0, 2, size=(n, 2, 1)).astype(np.int8)
    keys = p2.encode_subwindow_keys(hap, [(0, 1)], verify=False)
    Z, names, _ = p2.state_block(keys, np.arange(n))
    G = (hap[:, 0, 0] + hap[:, 1, 0]).astype(float)
    nz = [j for j in range(Z.shape[1]) if Z[:, j].std() > 0]
    ok = bool(nz) and all(np.allclose(np.corrcoef(Z[:, j], G)[0, 1] ** 2, 1.0) for j in nz)
    check("single_marker_state_is_function_of_G", ok,
          "ncols=%d nonconstant=%d" % (Z.shape[1], len(nz)))

def t_label_leakage():
    rng = np.random.RandomState(9)
    keys = np.array([["0" * 16, "1" * 16]] * 60, dtype=object)
    try:
        hstate.StateDictionary.fit(keys, phenotype=np.zeros(60))
        check("label_leakage_sentinel", False, "fit accepted a phenotype argument")
    except hstate.LabelLeakageError:
        check("label_leakage_sentinel", True)

def t_unseen_and_missing():
    a = "0" * 16
    b = "1" + "0" * 15
    train = np.array([[a, a]] * 60 + [[b, b]] * 60, dtype=object)
    d = hstate.StateDictionary.fit(train, min_distinct_people=50, max_states=32)
    unseen = np.array([["1" * 16, a]], dtype=object)
    Z, cols = d.transform(unseen)
    check("unseen_state_to_OTHER", Z[0, cols.index("OTHER")] == 1.0)
    missk = "." + "0" * 15
    Zm, colsm = d.transform(np.array([[missk, a]], dtype=object))
    check("missing_not_counted_as_reference", Zm[0, colsm.index("OTHER")] == 1.0)

def t_state_label_permutation():
    rng = np.random.RandomState(13)
    n, m = 300, 16
    hap = rng.randint(0, 2, size=(n, 2, m)).astype(np.int8)
    uo, codes = p2.encode_subwindow_keys(hap, [(0, m)], verify=False)[0]
    keys = uo[codes]
    d = hstate.StateDictionary.fit(keys, min_distinct_people=5, max_states=32)
    Z1, c1 = d.transform(keys)
    order = list(range(len(d.states)))[::-1]
    d2 = d.reorder(order)
    Z2, c2 = d2.transform(keys)
    m1 = dict((c, Z1[:, j]) for j, c in enumerate(c1))
    ok = all(np.array_equal(m1[c], Z2[:, j]) for j, c in enumerate(c2))
    check("state_label_permutation_invariant", ok)

def t_ld_prune_deterministic():
    rng = np.random.RandomState(17)
    n = 500
    x = rng.normal(size=n)
    X = np.column_stack([x, x + 1e-6 * rng.normal(size=n), rng.normal(size=n)]).astype(np.float32)
    R = p2.corr_matrix(X, np.arange(n))
    k1 = p2.greedy_ld_prune(R, 0.9)
    k2 = p2.greedy_ld_prune(R, 0.9)
    check("ld_prune_deterministic_drops_duplicate", k1 == k2 and k1 == [0, 2], "kept=%s" % k1)

def t_constant_column_dropped():
    U, C, X, y = _toy(seed=21)
    U[:, 1 + C.shape[1] + 2] = 3.0
    n = X.shape[0]
    ff = p2.FoldFit(U, C.shape[1], np.arange(n))
    fit = ff.solve(np.arange(X.shape[1]), 0.0, 0.01)
    check("zero_variance_column_dropped", 2 not in set(fit["keep"].tolist())
          and fit["keep"].size == X.shape[1] - 1, "kept=%d" % fit["keep"].size)

def main():
    t_ridge_matches_direct()
    t_enet_matches_sklearn()
    t_cis_trans()
    t_copy_swap_invariance()
    t_single_marker_window()
    t_label_leakage()
    t_unseen_and_missing()
    t_state_label_permutation()
    t_ld_prune_deterministic()
    t_constant_column_dropped()
    print("---")
    print("FAILED: %s" % (FAILS if FAILS else "none"))
    return 1 if FAILS else 0

if __name__ == "__main__":
    sys.exit(main())
