import glob
import json
import numpy as np

import att_common_oracle as C

STRATA = ['ALL', 'COMMON', 'LOW', 'RARE', 'VRARE']

def main():
    rows = []
    for f in sorted(glob.glob(C.ATT + '/oracle_chr*.json')):
        j = json.load(open(f))
        rows += j['per_domain']
    print('도메인 %d개' % len(rows))

    out = []
    for s in STRATA:
        d = [r[s] for r in rows if s in r and 'oof_z' in r[s]]
        if not d:
            continue
        n = np.array([x['n'] for x in d], float)
        fz = np.array([x.get('flat_z', np.nan) for x in d], float)
        fzn = np.array([x.get('flat_z_null', np.nan) for x in d], float)
        sz = np.array([x['single_z'] for x in d], float)
        szn = np.array([x['single_z_null_mean'] for x in d], float)
        oz = np.array([x['oof_z'] for x in d], float)
        ozn = np.array([x['oof_z_null_mean'] for x in d], float)
        op = np.array([x['oof_p'] for x in d], float)
        out.append(dict(
            stratum=s, n_dom=len(d), n_var_median=float(np.median(n)),
            flat_obs=float(np.nanmedian(fz)),
            flat_null=float(np.nanmedian(fzn)),
            single_obs=float(np.nanmedian(sz)),
            single_null=float(np.nanmedian(szn)),
            single_excess=float(np.nanmedian(sz - szn)),
            oof_obs=float(np.nanmedian(oz)),
            oof_null=float(np.nanmedian(ozn)),
            oof_excess=float(np.nanmedian(oz - ozn)),
            oof_minus_single=float(np.nanmedian((oz - ozn) - (sz - szn))),
            frac_oof_p_lt005=float(np.nanmean(op < 0.05)),
            frac_oof_p_lt01=float(np.nanmean(op < 0.1))))

    h = '%-7s %6s %7s %9s %9s %11s %9s %9s %11s %12s'
    print()
    print(h % ('층', 'n도메', '변이중앙', 'flat', 'single관측', 'single초과',
               'oof관측', 'oof귀무', 'oof초과', 'oof-single'))
    for r in out:
        print(h % (r['stratum'], r['n_dom'], '%.0f' % r['n_var_median'],
                   '%.3f' % r['flat_obs'], '%.3f' % r['single_obs'],
                   '%+.3f' % r['single_excess'], '%.3f' % r['oof_obs'],
                   '%.3f' % r['oof_null'], '%+.3f' % r['oof_excess'],
                   '%+.3f' % r['oof_minus_single']))
    print()
    print('%-7s %-22s %-22s' % ('층', 'oof p<0.05 도메인 비율', 'oof p<0.10 비율'))
    for r in out:
        print('%-7s %-22s %-22s' % (r['stratum'],
                                    '%.1f%%' % (100 * r['frac_oof_p_lt005']),
                                    '%.1f%%' % (100 * r['frac_oof_p_lt01'])))
    p = C.ATT + '/oracle_summary.json'
    json.dump(dict(n_domains=len(rows), nperm=C.NPERM, trait='tchl',
                   alphas=C.ALPHAS, result=out), open(p, 'w'),
              indent=1, ensure_ascii=False)
    print('\n->', p)

if __name__ == '__main__':
    main()
