import os
import sys
import json
import numpy as np
from sklearn.metrics import adjusted_rand_score as ari

import f3f6_common as C

N_F46 = int(os.environ.get('N_F46', '300'))

def main():
    os.makedirs(C.OUT, exist_ok=True)
    if not os.path.exists(C.CACHE + '/master.npz'):
        meta = C.build_cache()
        print('cache', json.dumps(meta))
    d = C.load()
    sel = C.pick_domains(d, N_F46, random_only=True)
    np.save(C.OUT + '/f46_domains.npy', sel)
    print('f46_domains', len(sel), 'forced_included',
          int(d['dom_forced'][sel].sum()))

    di = d['dom_idx']
    combos = ['spectral_eucl', 'ward_eucl', 'average_eucl', 'average_bio',
              'spectral_bio', 'dbscan_bio']
    rows = []
    for g in range(len(d['udom'])):
        rid = np.where(di == g)[0]
        cc = (d['has_ccre'][rid] > 0).astype(int)
        tf = (d['n_tf'][rid] > 0).astype(int)
        re2 = (d['n_re2g'][rid] > 0).astype(int)
        ccl = d['has_ccre'][rid].astype(int)
        for cb in combos:
            lk = d['mod_' + cb][rid]
            m = lk >= 0
            if m.sum() < 4 or len(np.unique(lk[m])) < 2:
                continue
            rows.append([int(g), int(d['dom_forced'][g]), cb, int(m.sum()),
                         int(len(np.unique(lk[m]))),
                         float(ari(lk[m], cc[m])), float(ari(lk[m], tf[m])),
                         float(ari(lk[m], re2[m])),
                         float(ari(lk[m], ccl[m])),
                         float(cc[m].mean()), float(tf[m].mean()),
                         float(re2[m].mean())])
    hdr = ['domain_idx', 'forced', 'combo', 'n_rows', 'k_modules',
           'ari_vs_has_ccre', 'ari_vs_has_tf', 'ari_vs_has_re2g',
           'ari_vs_ccre_count', 'frac_has_ccre', 'frac_has_tf',
           'frac_has_re2g']
    with open(C.OUT + '/f3b_label_vs_annotation.csv', 'w') as fh:
        fh.write(','.join(hdr) + '\n')
        for r in rows:
            fh.write(','.join([str(x) for x in r]) + '\n')
    R = np.array([[r[5], r[6], r[7]] for r in rows if r[2] == 'spectral_eucl'],
                 dtype=float)
    with open(C.OUT + '/prep.done', 'w') as fh:
        json.dump(dict(n_f46_domains=int(len(sel)),
                       n_rows_annot=len(rows),
                       spectral_eucl_median_ari_vs_has_ccre=float(
                           np.median(R[:, 0])) if len(R) else None,
                       spectral_eucl_median_ari_vs_has_tf=float(
                           np.median(R[:, 1])) if len(R) else None,
                       spectral_eucl_median_ari_vs_has_re2g=float(
                           np.median(R[:, 2])) if len(R) else None), fh,
                  indent=1)

if __name__ == '__main__':
    main()
