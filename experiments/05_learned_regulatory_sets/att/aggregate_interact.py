import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import glob
import json
import numpy as np

BASE = _config_path('${PROJECT_ROOT}/work/fset/att/interact')

def main():
    rows = []
    for f in sorted(glob.glob(BASE + '/int_chr*.json')):
        rows += json.load(open(f))['per_domain']
    print('도메인 %d개' % len(rows))
    if not rows:
        return
    nv = np.array([r['n_var'] for r in rows], float)
    npr = np.array([r['n_pairs'] for r in rows], float)
    print('단위당 변이 중앙 %.0f (쌍 %.0f개)' % (np.median(nv), np.median(npr)))
    print()
    res = []
    for nm in ('MAIN', 'MAIN_INT'):
        d = [r[nm] for r in rows if nm in r]
        o = np.array([x['obs'] for x in d])
        n = np.array([x['null_mean'] for x in d])
        e = np.array([x['excess'] for x in d])
        p = np.array([x['p'] for x in d])
        res.append(dict(arm=nm, n=len(d), obs=float(np.nanmedian(o)),
                        null=float(np.nanmedian(n)),
                        excess=float(np.nanmedian(e)),
                        excess_p75=float(np.nanpercentile(e, 75)),
                        frac_p005=float(np.nanmean(p < 0.05)),
                        frac_p01=float(np.nanmean(p < 0.1))))
    h = '%-10s %6s %9s %9s %10s %12s %11s %10s'
    print(h % ('팔', 'n', '관측중앙', '귀무중앙', '초과중앙', '초과p75',
               'p<0.05', 'p<0.10'))
    for r in res:
        print(h % (r['arm'], r['n'], '%.3f' % r['obs'], '%.3f' % r['null'],
                   '%+.4f' % r['excess'], '%+.4f' % r['excess_p75'],
                   '%.1f%%' % (100 * r['frac_p005']),
                   '%.1f%%' % (100 * r['frac_p01'])))
    pair = [(r['MAIN']['excess'], r['MAIN_INT']['excess'])
            for r in rows if 'MAIN' in r and 'MAIN_INT' in r]
    if pair:
        a = np.array([x[0] for x in pair])
        b = np.array([x[1] for x in pair])
        dlt = b - a
        print()
        print('도메인별 (MAIN_INT 초과 - MAIN 초과):')
        print('  중앙 %+.4f, 평균 %+.4f, >0 비율 %.1f%%, p90 %+.4f'
              % (np.nanmedian(dlt), np.nanmean(dlt), 100 * np.nanmean(dlt > 0),
                 np.nanpercentile(dlt, 90)))
    json.dump(dict(n_domains=len(rows), result=res),
              open(BASE + '/interact_summary.json', 'w'), indent=1,
              ensure_ascii=False)
    print('\n->', BASE + '/interact_summary.json')

if __name__ == '__main__':
    main()
