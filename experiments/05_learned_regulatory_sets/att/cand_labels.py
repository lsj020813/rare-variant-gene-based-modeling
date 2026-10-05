import json
import numpy as np

import att_common as C

OUT = C.ATT + '/cand_labels.npz'

def collapse2(lab, di, ndom):
    v = np.full(len(di), -1, dtype=np.int8)
    for g in range(ndom):
        rid = np.where(di == g)[0]
        l = lab[rid]
        m = l >= 0
        if m.sum() < 2:
            continue
        u, cnt = np.unique(l[m], return_counts=True)
        if len(u) < 2:
            continue
        mj = u[int(np.argmax(cnt))]
        w = np.where(l == mj, 0, 1).astype(np.int8)
        w[~m] = -1
        v[rid] = w
    return v

def median_split(x, di, ndom, rng=None):
    v = np.full(len(di), -1, dtype=np.int8)
    for g in range(ndom):
        rid = np.where(di == g)[0]
        a = x[rid]
        m = np.isfinite(a)
        if m.sum() < 4:
            continue
        t = np.median(a[m])
        w = np.where(a > t, 0, 1).astype(np.int8)
        if (w[m] == 0).all() or (w[m] == 1).all():
            continue
        w[~m] = -1
        v[rid] = w
    return v

def main():
    d = C.master()
    di = np.asarray(d['dom_idx'])
    X = np.asarray(d['X'], dtype=np.float32)
    xc = list(np.asarray(d['X_cols']))
    maf = np.asarray(d['maf'], dtype=np.float64)
    r2 = np.asarray(d['r2'], dtype=np.float64)
    hcc = np.asarray(d['has_ccre'], dtype=np.float64)
    ntf = np.asarray(d['n_tf'], dtype=np.float64)
    dcen = np.asarray(d['dist_center'], dtype=np.float64)
    nd = C.NDOM
    rng = np.random.default_rng(C.SEED + 909)

    def col(name):
        return X[:, xc.index(name)].astype(np.float64)

    parts, names, kinds = [], [], []

    def add(nm, v, kind):
        parts.append(v)
        names.append(nm)
        kinds.append(kind)
        ok = np.array([(v[di == g] >= 0).any() for g in range(nd)])
        print('  %-18s %-14s 유효도메인 %4d  모듈0 비율 %.3f'
              % (nm, kind, ok.sum(),
                 float((v[v >= 0] == 0).mean()) if (v >= 0).any() else 0.0),
              flush=True)

    for k in ('mod_spectral_eucl', 'mod_ward_eucl', 'mod_average_eucl',
              'mod_average_bio', 'mod_spectral_bio', 'mod_dbscan_bio'):
        add(k.replace('mod_', ''), collapse2(np.asarray(d[k]), di, nd),
            'annot-cluster')

    add('maf', median_split(np.log10(np.clip(maf, 1e-8, None)), di, nd),
        'freq')
    add('r2', median_split(r2, di, nd), 'imputation')
    add('dist_center', median_split(dcen, di, nd), 'position')
    add('gpn', median_split(col('gpn_score_rz'), di, nd), 'seq-model')
    add('tf_n', median_split(col('tf_n_rz'), di, nd), 'annot-density')
    add('ccre_any', median_split(hcc + 1e-6 * rng.standard_normal(len(hcc)),
                                 di, nd), 'annot-binary')
    add('random', median_split(rng.standard_normal(len(di)), di, nd),
        'negative-control')

    L = np.stack(parts).astype(np.int8)
    ok = np.array([(L[:, di == g] >= 0).all(1).any() for g in range(nd)])
    okall = np.ones(nd, dtype=bool)
    for p in range(len(parts)):
        v = parts[p]
        for g in range(nd):
            if not (v[di == g] >= 0).any():
                okall[g] = False
    print('\n분할 %d개, 전 분할 공통 유효 도메인 %d' % (len(parts), okall.sum()))
    np.savez_compressed(OUT, labels=L, names=np.array(names),
                        kinds=np.array(kinds), dom_ok=okall)
    json.dump(dict(n_partitions=len(parts), names=names, kinds=kinds,
                   n_dom_ok=int(okall.sum()), seed=C.SEED),
              open(C.ATT + '/cand_labels_meta.json', 'w'),
              indent=1, ensure_ascii=False)
    print('->', OUT)

if __name__ == '__main__':
    main()
