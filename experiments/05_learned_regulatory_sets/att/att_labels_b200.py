import json
import numpy as np
import att_common_b200 as C

def two_level(lab, di):
    n = len(di)
    lab2 = np.full(n, -1, dtype=np.int8)
    ok = np.zeros(C.NDOM, dtype=bool)
    nk = np.zeros((C.NDOM, 2), dtype=np.int32)
    merged, single, kvals = [], [], []
    for g in range(C.NDOM):
        rid = np.where(di == g)[0]
        l = lab[rid]
        m = l >= 0
        u, cnt = np.unique(l[m], return_counts=True)
        kvals.append(int(len(u)))
        if len(u) < 2:
            single.append(int(g))
            continue
        if len(u) > 2:
            merged.append(int(g))
        mj = u[int(np.argmax(cnt))]
        v = np.where(l == mj, 0, 1).astype(np.int8)
        v[~m] = -1
        lab2[rid] = v
        ok[g] = True
        nk[g, 0] = int((v == 0).sum())
        nk[g, 1] = int((v == 1).sum())
    return lab2, ok, nk, merged, single, kvals

def shuffle_labels(lab2, di, ok, seed):
    rng = np.random.default_rng(seed)
    out = lab2.copy()
    for g in range(C.NDOM):
        if not ok[g]:
            continue
        rid = np.where(di == g)[0]
        sel = rid[lab2[rid] >= 0]
        v = lab2[sel].copy()
        rng.shuffle(v)
        out[sel] = v
    return out

def features(lab2, d, di, ok):
    maf = np.asarray(d['maf'], dtype=np.float64)
    r2 = np.asarray(d['r2'], dtype=np.float64)
    nre = np.asarray(d['n_re2g'], dtype=np.float64)
    hcc = np.asarray(d['has_ccre'], dtype=np.float64)
    X = np.asarray(d['X'], dtype=np.float32)
    raw = np.full((C.NDOM, 2, 6), np.nan, dtype=np.float32)
    cm = np.full((C.NDOM, 2, X.shape[1]), np.nan, dtype=np.float32)
    for g in range(C.NDOM):
        if not ok[g]:
            continue
        rid = np.where(di == g)[0]
        lv = lab2[rid]
        ntot = float((lv >= 0).sum())
        for k in (0, 1):
            s = rid[lv == k]
            if len(s) == 0:
                continue
            mk = np.clip(maf[s], 1e-6, None)
            raw[g, k, 0] = np.log10(len(s))
            raw[g, k, 1] = len(s) / ntot
            raw[g, k, 2] = float(np.nanmean(np.log10(mk)))
            raw[g, k, 3] = float(np.nanmean(r2[s]))
            raw[g, k, 4] = float(np.nanmean(nre[s] > 0))
            raw[g, k, 5] = float(np.nanmean(hcc[s] > 0))
            cm[g, k] = X[s].mean(axis=0)
    return raw, cm

def main():
    C.ensure_dirs()
    d = C.master()
    di = np.asarray(d['dom_idx'])
    lab = np.asarray(d['mod_spectral_eucl'])
    lab2, ok, nk, merged, single, kvals = two_level(lab, di)
    L = np.zeros((C.NSHUF + 1, len(di)), dtype=np.int8)
    L[0] = lab2
    for i in range(1, C.NSHUF + 1):
        L[i] = shuffle_labels(lab2, di, ok, C.SEED + i)
    raws, cms = [], []
    for i in range(C.NSHUF + 1):
        a, b = features(L[i], d, di, ok)
        raws.append(a)
        cms.append(b)
    np.savez_compressed(
        C.ATT + '/att_labels.npz', labels=L, dom_ok=ok, nk=nk,
        raw=np.stack(raws), cmean=np.stack(cms),
        x_cols=np.array(['log10_nvar', 'size_frac', 'mean_log10_maf',
                         'mean_r2', 'frac_re2g', 'frac_ccre_diag']),
        cset_cols=np.asarray(d['X_cols']))
    meta = dict(n_dom_usable=int(ok.sum()), n_dom_single_cluster=len(single),
                n_dom_merged_from_K_gt2=len(merged), merged_domains=merged,
                k_hist={str(k): int((np.array(kvals) == k).sum())
                        for k in sorted(set(kvals))},
                n_shuffles=C.NSHUF, seed=C.SEED,
                nk_major_median=float(np.median(nk[ok, 0])),
                nk_minor_median=float(np.median(nk[ok, 1])))
    with open(C.ATT + '/att_labels_meta.json', 'w') as f:
        json.dump(meta, f, indent=1)
    print(json.dumps(meta, indent=1)[:1500])

if __name__ == '__main__':
    main()
