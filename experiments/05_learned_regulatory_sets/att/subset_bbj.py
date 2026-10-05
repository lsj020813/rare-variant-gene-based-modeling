import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import sys
import csv
import json
import numpy as np
from scipy.stats import norm

sys.path.insert(0, _config_path('${PROJECT_ROOT}/work/fset/att'))
import att_common as C

FITD = C.ATT + '/' + (sys.argv[1] if len(sys.argv) > 1 else 'lamfix')
BBJ = C.ATT + '/att_bbj_tc_domains.csv'

def bh(p):
    p = np.asarray(p, dtype=np.float64)
    m = np.isfinite(p)
    q = np.full(len(p), np.nan)
    pp = p[m]
    o = np.argsort(pp)
    n = len(pp)
    r = np.empty(n)
    r[o] = np.arange(1, n + 1)
    qq = pp * n / r
    qs = np.minimum.accumulate(qq[o][::-1])[::-1]
    out = np.empty(n)
    out[o] = qs
    q[m] = np.minimum(out, 1.0)
    return q

def summarize(zatt, zflat, tag):
    q = bh(2.0 * norm.sf(np.abs(zatt)))
    return dict(tag=tag, n=len(zatt),
                n_sig=int(np.nansum(q < 0.05)),
                mean_delta=float(np.nanmean(np.abs(zatt) - np.abs(zflat))),
                max_abs_z=float(np.nanmax(np.abs(zatt))))

def pval(obs, null):
    null = np.asarray(null, dtype=float)
    return (1 + int((null >= obs).sum())) / (len(null) + 1)

def main():
    g = list(csv.DictReader(open(FITD + '/att_gene_results.csv')))
    zflat = np.array([float(r['z_FLAT']) for r in g])
    zatt = np.array([float(r['z_ATT']) for r in g])
    zint = np.array([float(r['z_INT']) for r in g])

    bb = {int(r['domain_idx']): int(r['bbj_tc_domain'])
          for r in csv.DictReader(open(BBJ))}
    sel = np.array([bb.get(i, 0) == 1 for i in range(len(g))])
    print('BBJ TC 플래그 도메인 %d / %d' % (sel.sum(), len(g)))

    nulls = sorted(__import__('glob').glob(FITD + '/private/null_z_*.npy'))
    ndi = sorted(__import__('glob').glob(FITD
                                         + '/private/null_delta_int_*.npy'))
    print('셔플 %d회 (ATT), %d회 (INT)' % (len(nulls), len(ndi)))

    res = {}
    for scope, m in (('전체 1000', np.ones(len(g), bool)), ('BBJ 22', sel)):
        obs_a = summarize(zatt[m], zflat[m], 'ATT')
        obs_i = summarize(zint[m], zflat[m], 'INT')
        na = [summarize(np.load(f)[m], zflat[m], 'ATT') for f in nulls]
        ni = [summarize(np.load(f)[m] + np.abs(zflat[m]), zflat[m], 'INT')
              for f in ndi]
        row = {}
        for arm, obs, nl in (('ATT', obs_a, na), ('INT', obs_i, ni)):
            for k in ('n_sig', 'mean_delta', 'max_abs_z'):
                nv = [x[k] for x in nl]
                row['%s_%s' % (arm, k)] = dict(
                    obs=obs[k], null_mean=float(np.mean(nv)),
                    null_min=float(np.min(nv)), null_max=float(np.max(nv)),
                    p=pval(obs[k], nv))
        res[scope] = dict(n_genes=int(m.sum()), stats=row)

    print()
    hdr = '%-10s %-22s %10s %10s %16s %8s'
    print(hdr % ('범위', '통계량', '관측', '귀무평균', '귀무범위', '순열 p'))
    for scope, d in res.items():
        for k, v in d['stats'].items():
            print(hdr % (scope, k, '%.4f' % v['obs'], '%.4f' % v['null_mean'],
                         '[%.3f, %.3f]' % (v['null_min'], v['null_max']),
                         '%.4f' % v['p']))
        print()

    with open(FITD + '/subset_bbj.json', 'w') as fh:
        json.dump(dict(fit_dir=FITD, n_bbj=int(sel.sum()),
                       n_shuffles=len(nulls),
                       p_floor=1.0 / (len(nulls) + 1),
                       fdr_within_subset=True, result=res),
                  fh, indent=1, ensure_ascii=False)
    print('-> %s/subset_bbj.json' % FITD)

if __name__ == '__main__':
    main()
