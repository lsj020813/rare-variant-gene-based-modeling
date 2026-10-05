import os as _cfg_os
import math as _cfg_math

def _config_number(name, cast=int, positive=False):
    value = cast(_cfg_os.environ[name])
    if not _cfg_math.isfinite(value) or (positive and value <= 0):
        raise ValueError(name + " must be finite" + (" and positive" if positive else ""))
    return value

import os
import glob
import json
import numpy as np

import f3f6_common as C

KTOP = 21

def main():
    d = C.load()
    rows = []
    hdr = None
    for p in sorted(glob.glob(C.OUT + '/f6_raw_chr*.csv')):
        with open(p) as fh:
            h = fh.readline().rstrip('\n').split(',')
            if hdr is None:
                hdr = h
            for line in fh:
                q = line.rstrip('\n').split(',')
                if len(q) == len(h):
                    rows.append(q)
    I = dict((c, i) for i, c in enumerate(hdr))
    n = len(rows)
    car = np.array([float(r[I['n_carriers']]) for r in rows])
    bvar = np.array([float(r[I['burden_var']]) for r in rows])
    grade = np.where((car >= 100) & (bvar > 0), 'USABLE',
                     np.where(car >= 10, 'MARGINAL', 'UNUSABLE'))
    with open(C.OUT + '/module_statistical_usability.csv', 'w') as fh:
        fh.write(','.join(hdr[:16] + ['usability_grade']) + '\n')
        for i, r in enumerate(rows):
            fh.write(','.join(r[:16] + [grade[i]]) + '\n')

    K = np.array([[float(r[I['K%d' % j]]) for j in range(KTOP)] for r in rows])
    forced = np.array([int(r[I['forced']]) for r in rows])
    out = []
    for scope, msk in [('all_sampled', np.ones(n, dtype=bool)),
                       ('rand_only', forced == 0),
                       ('forced_only', forced == 1)]:
        if msk.sum() == 0:
            continue
        cs = car[msk]
        ks = K[msk]
        ext = cs * _config_number("N_SAMPLES", float, True) / 5000.0
        row = dict(
            scope=scope, n_modules=int(msk.sum()),
            n_domains=int(len(set([rows[i][I['domain_idx']]
                                   for i in np.where(msk)[0]]))),
            median_n_var=float(np.median([float(rows[i][I['n_var']])
                                          for i in np.where(msk)[0]])),
            median_carriers_sub5000=float(np.median(cs)),
            median_carriers_extrap=float(np.median(ext)),
            n_ge5=int((cs >= 5).sum()), n_ge10=int((cs >= 10).sum()),
            n_ge50=int((cs >= 50).sum()), n_ge100=int((cs >= 100).sum()),
            frac_ge100=float((cs >= 100).mean()),
            n_usable=int((grade[msk] == 'USABLE').sum()),
            n_marginal=int((grade[msk] == 'MARGINAL').sum()),
            n_unusable=int((grade[msk] == 'UNUSABLE').sum()),
            median_burden_var=float(np.median(bvar[msk])),
            median_cocarry_jaccard=float(np.nanmedian(
                [float(rows[i][I['cocarry_jaccard']] or 'nan')
                 for i in np.where(msk)[0]])),
            median_mean_K_per_carrier=float(np.median(
                [float(rows[i][I['mean_K_per_carrier']])
                 for i in np.where(msk)[0]])),
            frac_carriers_with_K_ge2=float(
                ks[:, 2:].sum() / max(ks[:, 1:].sum(), 1.0)),
            frac_carriers_with_K_ge3=float(
                ks[:, 3:].sum() / max(ks[:, 1:].sum(), 1.0)))
        out.append(row)
    keys = list(out[0].keys())
    with open(C.OUT + '/module_context_recurrence.csv', 'w') as fh:
        fh.write(','.join(keys) + '\n')
        for r in out:
            fh.write(','.join([str(r[k]) for k in keys]) + '\n')
        fh.write('\n# per-module K distribution (recurrence of the module '
                 'context within a person), pooled over modules\n')
        fh.write('scope,' + ','.join(['K%d_persons' % j
                                      for j in range(KTOP)]) + '\n')
        for scope, msk in [('all_sampled', np.ones(n, dtype=bool)),
                           ('rand_only', forced == 0),
                           ('forced_only', forced == 1)]:
            if msk.sum() == 0:
                continue
            fh.write(scope + ',' + ','.join(
                [str(int(K[msk][:, j].sum())) for j in range(KTOP)]) + '\n')
    with open(C.OUT + '/f6.done', 'w') as fh:
        json.dump(dict(n_modules=n, grades=dict(
            USABLE=int((grade == 'USABLE').sum()),
            MARGINAL=int((grade == 'MARGINAL').sum()),
            UNUSABLE=int((grade == 'UNUSABLE').sum()))), fh, indent=1)

if __name__ == '__main__':
    main()
