import os
import json
import numpy as np
from sklearn.cluster import SpectralClustering
from sklearn.metrics import adjusted_rand_score as ari

import f3f6_common as C

GENO = C.OUT + '/geno'
BOOT = C.BOOT
BOOT_ORD = 200
NSHUF = 20
RULE = ('F4 PASS <=> (global_pattern in {B,C}) AND '
        '(median ARI(residual, original) >= 0.4) AND '
        '(residual F1 tendency kept: >=50% of domains with Hopkins_obs > all '
        '20 shuffles AND 1nn_obs < all 20 shuffles); '
        'global_pattern: A if frac_module_A >= 0.5, B if frac_module_A <= 0.10,'
        ' else C; module A if obs mean r2 > null q97.5')

def hopkins(Xe, rng, frac=0.1):
    n, p = Xe.shape
    m = max(5, int(np.ceil(frac * n)))
    m = min(m, n - 1)
    lo = Xe.min(0)
    hi = Xe.max(0)
    ridx = rng.choice(n, m, replace=False)
    U = rng.uniform(lo, hi, size=(m, p))
    d2 = ((Xe[:, None, :] - U[None, :, :]) ** 2).sum(-1)
    u = np.sqrt(d2.min(0))
    d2r = ((Xe[:, None, :] - Xe[ridx][None, :, :]) ** 2).sum(-1)
    d2r[ridx, np.arange(m)] = np.inf
    w = np.sqrt(d2r.min(0))
    s = u.sum() + w.sum()
    return float(u.sum() / s) if s > 0 else np.nan

def nn1_median(Xe):
    d2 = ((Xe[:, None, :] - Xe[None, :, :]) ** 2).sum(-1)
    np.fill_diagonal(d2, np.inf)
    return float(np.median(np.sqrt(d2.min(1))))

def cmds(Dm, k):
    n = Dm.shape[0]
    J = np.eye(n) - 1.0 / n
    B = -0.5 * J @ (Dm ** 2) @ J
    w, V = np.linalg.eigh(B)
    o = np.argsort(w)[::-1][:k]
    w = np.clip(w[o], 0, None)
    return V[:, o] * np.sqrt(w)[None, :]

def main():
    d = C.load()
    lab = d['mod_spectral_eucl']
    di = d['dom_idx']
    sel = np.load(C.OUT + '/f46_domains.npy')
    X = d['X']
    logmaf = np.log10(np.maximum(d['maf'], 1e-6))

    mod_rows = []
    res_rows = []
    for g in sel:
        p = GENO + '/dom%06d.npz' % g
        if not os.path.exists(p):
            continue
        z = np.load(p)
        keep = z['local_rows']
        R2 = np.asarray(z['r2'], dtype=np.float64)
        ok = np.asarray(z['nonzero_var'], dtype=bool)
        rid = np.where(di == g)[0][keep]
        n = len(rid)
        if n < 6:
            continue
        R2f = np.nan_to_num(R2, nan=0.0)
        R2f[~ok, :] = 0.0
        R2f[:, ~ok] = 0.0
        np.fill_diagonal(R2f, ok.astype(np.float64))
        okf = ok.astype(np.float64)
        with np.errstate(invalid='ignore'):
            ldctx = (R2f.sum(1) - okf) / np.maximum(okf.sum() - 1.0, 1.0)
        ldctx = np.where(ok, ldctx, np.nan)
        strat = C.strata_for(d, rid, 'LD', ld_ctx=ldctx)
        lk = lab[rid]

        for u in np.unique(lk):
            if u < 0:
                continue
            mm = lk == u
            if int((mm & ok).sum()) < 2:
                continue
            zo = np.zeros((1, n))
            zo[0, mm] = 1.0
            idx = C.matched_null(strat, mm, boot=BOOT,
                                 seed=C.det_seed('f4', int(g), int(u)))
            Z = np.zeros((BOOT, n))
            Z[np.arange(BOOT)[:, None], idx] = 1.0

            def quad(Zm):
                M = Zm @ okf
                S = ((Zm * okf) @ R2f * (Zm * okf)).sum(1)
                with np.errstate(invalid='ignore', divide='ignore'):
                    mr = (S - M) / (M * (M - 1.0))
                    varl = S / M - 1.0
                    meff = 1.0 + (M - 1.0) * (1.0 - varl / M)
                return M, mr, np.where(M > 0, meff / M, np.nan)

            Mo, mro, mfo = quad(zo)
            Mn, mrn, mfn = quad(Z)
            gm = np.where(mm & ok)[0]
            iu = np.triu_indices(len(gm), 1)
            vo = R2[np.ix_(gm, gm)][iu]
            vo = vo[np.isfinite(vo)]
            medo = float(np.median(vo)) if len(vo) else np.nan
            maxo = float(np.max(vo)) if len(vo) else np.nan
            medn = np.full(BOOT_ORD, np.nan)
            maxn = np.full(BOOT_ORD, np.nan)
            for b in range(BOOT_ORD):
                s = idx[b]
                s = s[ok[s]]
                if len(s) < 2:
                    continue
                iu2 = np.triu_indices(len(s), 1)
                vv = R2[np.ix_(s, s)][iu2]
                vv = vv[np.isfinite(vv)]
                if len(vv):
                    medn[b] = np.median(vv)
                    maxn[b] = np.max(vv)
            fin = np.isfinite(mrn)
            q975 = float(np.quantile(mrn[fin], 0.975)) if fin.any() else np.nan
            q025 = float(np.quantile(mrn[fin], 0.025)) if fin.any() else np.nan
            pat = 'A' if (np.isfinite(mro[0]) and np.isfinite(q975)
                          and mro[0] > q975) else 'B'
            mod_rows.append([int(g), int(d['dom_forced'][g]), int(u),
                             int(Mo[0]), float(mro[0]),
                             float(np.nanmean(mrn)), q025, q975,
                             int(np.nansum(mrn >= mro[0])), medo,
                             float(np.nanmean(medn)), maxo,
                             float(np.nanmean(maxn)), float(mfo[0]),
                             float(np.nanmean(mfn)),
                             int(np.nansum(mfn <= mfo[0])), pat,
                             float((Z[:, mm].sum(1) / max(int(mm.sum()), 1))
                                   .mean())])

        iu = np.triu_indices(n, 1)
        Dx = np.sqrt(np.maximum(((X[rid][:, None, :] - X[rid][None, :, :])
                                 ** 2).sum(-1), 0.0)) \
            if n <= 420 else None
        if Dx is None:
            continue
        y = -Dx[iu]
        x1 = np.log10(1.0 + np.abs(d['pos38'][rid][:, None]
                                   - d['pos38'][rid][None, :]))[iu]
        x2 = np.nan_to_num(R2[iu], nan=0.0)
        x3 = np.abs(logmaf[rid][:, None] - logmaf[rid][None, :])[iu]
        A = np.column_stack([np.ones_like(y), x1, x2, x3])
        beta, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
        pred = A @ beta
        e = y - pred
        ss = float(((y - y.mean()) ** 2).sum())
        r2reg = 1.0 - float((e ** 2).sum()) / ss if ss > 0 else np.nan
        E = np.zeros((n, n))
        E[iu] = e
        E = E + E.T
        Dres = float(e.max()) - E
        np.fill_diagonal(Dres, 0.0)
        Dres = np.maximum(Dres, 0.0)
        k = int(min(10, n - 1))
        try:
            sc = SpectralClustering(n_clusters=2,
                                    affinity='precomputed_nearest_neighbors',
                                    n_neighbors=k, assign_labels='kmeans',
                                    random_state=C.SEED, n_init=10)
            lres = sc.fit_predict(Dres)
        except Exception:
            lres = np.zeros(n, dtype=int)
        cc = (d['has_ccre'][rid] > 0).astype(int)
        tf = (d['n_tf'][rid] > 0).astype(int)
        re2 = (d['n_re2g'][rid] > 0).astype(int)
        emb = cmds(Dres, int(min(5, n - 1)))
        rng = np.random.default_rng(C.det_seed('f4res', int(g)))
        h_o = hopkins(emb, rng)
        n_o = nn1_median(emb)
        hs = []
        ns = []
        for b in range(NSHUF):
            Es = emb.copy()
            for j in range(Es.shape[1]):
                Es[:, j] = Es[rng.permutation(n), j]
            hs.append(hopkins(Es, rng))
            ns.append(nn1_median(Es))
        hs = np.array(hs)
        ns = np.array(ns)
        res_rows.append([int(g), int(d['dom_forced'][g]), n,
                         float(ari(lk, lres)), float(ari(lk, cc)),
                         float(ari(lres, cc)), float(ari(lk, tf)),
                         float(ari(lk, re2)), r2reg, float(beta[1]),
                         float(beta[2]), float(beta[3]), h_o,
                         float(hs.min()), float(hs.max()),
                         int((hs < h_o).sum()), n_o, float(ns.min()),
                         float(ns.max()), int((ns > n_o).sum()),
                         int(len(np.unique(lk))), int(len(np.unique(lres)))])

    hdr = ['domain_idx', 'forced', 'module', 'n_var_ok', 'obs_mean_r2',
           'null_mean_r2', 'null_q025', 'null_q975', 'n_null_ge_obs',
           'obs_median_r2', 'null_median_r2_mean200', 'obs_max_r2',
           'null_max_r2_mean200', 'obs_Meff_over_M', 'null_Meff_over_M_mean',
           'n_null_Meff_le_obs', 'pattern', 'mean_null_module_overlap']
    with open(C.OUT + '/ld_confounding.csv', 'w') as fh:
        fh.write(','.join(hdr) + '\n')
        for r in mod_rows:
            fh.write(','.join([str(x) for x in r]) + '\n')
    hdr2 = ['domain_idx', 'forced', 'n_var', 'ari_residual_vs_original',
            'ari_original_vs_has_ccre', 'ari_residual_vs_has_ccre',
            'ari_original_vs_has_tf', 'ari_original_vs_has_re2g',
            'reg_r2_sim_on_dist_ld_maf', 'beta_logdist', 'beta_r2',
            'beta_dmaf', 'resid_hopkins', 'resid_hopkins_shuf_min',
            'resid_hopkins_shuf_max', 'n_shuf_below_obs', 'resid_nn1',
            'resid_nn1_shuf_min', 'resid_nn1_shuf_max', 'n_shuf_above_obs',
            'k_original', 'k_residual']
    with open(C.OUT + '/residual_regulatory_similarity.csv', 'w') as fh:
        fh.write(','.join(hdr2) + '\n')
        for r in res_rows:
            fh.write(','.join([str(x) for x in r]) + '\n')

    M = np.array([r[16] == 'A' for r in mod_rows], dtype=bool)
    fa = float(M.mean()) if len(M) else float('nan')
    gp = 'A' if fa >= 0.5 else ('B' if fa <= 0.10 else 'C')
    Rr = np.array([[r[3], r[15], r[19]] for r in res_rows], dtype=float) \
        if res_rows else np.zeros((0, 3))
    med_ari = float(np.nanmedian(Rr[:, 0])) if len(Rr) else float('nan')
    keep_f1 = float(np.mean((Rr[:, 1] >= NSHUF) & (Rr[:, 2] >= NSHUF))) \
        if len(Rr) else float('nan')
    verdict = 'PASS' if (gp in ('B', 'C') and med_ari >= 0.4
                         and keep_f1 >= 0.5) else 'FAIL'
    with open(C.OUT + '/f4_judgement.json', 'w') as fh:
        json.dump(dict(rule=RULE, n_modules=len(mod_rows),
                       n_domains=len(res_rows), frac_module_pattern_A=fa,
                       global_pattern=gp, median_ari_residual=med_ari,
                       frac_domains_residual_F1_kept=keep_f1,
                       boot=BOOT, boot_order=BOOT_ORD, verdict=verdict),
                  fh, indent=1)
    with open(C.OUT + '/f4.done', 'w') as fh:
        fh.write(verdict + '\n')

if __name__ == '__main__':
    main()
