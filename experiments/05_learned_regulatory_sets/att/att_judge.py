import os as _cfg_os
import math as _cfg_math

def _config_number(name, cast=int, positive=False):
    value = cast(_cfg_os.environ[name])
    if not _cfg_math.isfinite(value) or (positive and value <= 0):
        raise ValueError(name + " must be finite" + (" and positive" if positive else ""))
    return value

import os as _os
MDE_REFERENCE = _config_number("MDE_REFERENCE", float, True)
import csv
import json
import hashlib
import numpy as np
from scipy.stats import norm, binomtest
import att_common as C

RULES = {
 'primary_metric': 'Δ_g = |z_ATT| − |z_FLAT| (개인 out-of-fold)',
 'INT_appendix': ('A5 부록(상위 세션 지시, 표현형 읽기 전 접수): INT = FLAT + '
                  'b_major*b_minor + b_major^2 + b_minor^2 유전자별 회귀. '
                  'Δ_INT = |z_INT| − |z_FLAT| 에 ATT 와 동일한 3조건 적용 → '
                  'MODULE_INTERACTION_GAIN / NO_GAIN / NOT_ASSESSABLE'),
 'MODULE_ATTENTION_GAIN': ('평균 Δ > 라벨 셔플 20회 pooled 95분위 AND '
                           '짝지은 부호검정(Δ>0, 단측 α 0.025) AND '
                           'n_ATT > n_FLAT'),
 'NO_GAIN': '위 미충족이고 신호 존재(n_FLAT ≥ 1 또는 max|z_FLAT| ≥ 4)',
 'NOT_ASSESSABLE': 'n_FLAT = n_FIXED = n_ATT = 0 이고 max|z| < 4 전체',
 'null_primary': '셔플 20회의 평균 Δ 20개 값의 95분위 (선형보간)',
 'null_secondary': '셔플 20 × 1,000 유전자 Δ 20,000개 pooled 의 95분위',
 'fdr': 'BH 5%, 1,000 도메인, 팔별 독립 적용',
 'ref_MDE': MDE_REFERENCE,
 'ref_B_crossfit_z': [-4.22, -6.01],
}
ALPHA = 0.025
FDR = 0.05

def main():
    R = np.load(C.PRIV + '/att_raw.npz')
    z = dict((k[2:], R[k]) for k in R.files if k.startswith('z_'))
    q = {}
    for k, v in z.items():
        p = 2.0 * norm.sf(np.abs(v))
        q[k] = _bh(p)
    nsig = dict((k, int(np.nansum(v < FDR))) for k, v in q.items())
    dz = np.asarray(R['dz_obs'])
    good = np.isfinite(dz)
    mean_delta = float(np.nanmean(dz))
    nl = list(csv.DictReader(open(C.ATT + '/att_label_shuffle_null.csv')))
    null_means = np.array([float(r['mean_delta']) for r in nl])
    pooled = np.concatenate([np.load(C.PRIV + '/null_delta_%02d.npy'
                                     % int(r['shuffle'])) for r in nl])
    pooled = pooled[np.isfinite(pooled)]
    n95 = float(np.quantile(null_means, 0.95))
    p95 = float(np.quantile(pooled, 0.95))
    npos = int(np.nansum(dz[good] > 0))
    nnz = int(np.nansum(dz[good] != 0))
    st = binomtest(npos, nnz, 0.5, alternative='greater')
    sign_p = float(st.pvalue)
    c1 = mean_delta > n95
    c2 = sign_p < ALPHA
    c3 = nsig['ATT'] > nsig['FLAT']
    maxflat = float(np.nanmax(np.abs(z['FLAT'])))
    maxall = float(max(np.nanmax(np.abs(z[k])) for k in
                       ('FLAT', 'FIXED', 'ATT')))
    if nsig['FLAT'] == 0 and nsig['FIXED'] == 0 and nsig['ATT'] == 0 \
            and maxall < 4:
        verdict = 'NOT_ASSESSABLE'
    elif c1 and c2 and c3:
        verdict = 'MODULE_ATTENTION_GAIN'
    elif nsig['FLAT'] >= 1 or maxflat >= 4:
        verdict = 'NO_GAIN'
    else:
        verdict = 'NOT_ASSESSABLE'
    dzi = np.asarray(R['dz_int_obs'])
    gi = np.isfinite(dzi)
    mean_delta_int = float(np.nanmean(dzi))
    null_means_int = np.array([float(r['mean_delta_int']) for r in nl])
    pooled_int = np.concatenate([np.load(C.PRIV + '/null_delta_int_%02d.npy'
                                         % int(r['shuffle'])) for r in nl])
    pooled_int = pooled_int[np.isfinite(pooled_int)]
    n95i = float(np.quantile(null_means_int, 0.95))
    p95i = float(np.quantile(pooled_int, 0.95))
    nposi = int(np.nansum(dzi[gi] > 0))
    nnzi = int(np.nansum(dzi[gi] != 0))
    sign_pi = float(binomtest(nposi, nnzi, 0.5, alternative='greater').pvalue)
    i1 = mean_delta_int > n95i
    i2 = sign_pi < ALPHA
    i3 = nsig['INT'] > nsig['FLAT']
    if nsig['FLAT'] == 0 and nsig['FIXED'] == 0 and nsig['INT'] == 0 \
            and maxall < 4:
        verdict_int = 'NOT_ASSESSABLE'
    elif i1 and i2 and i3:
        verdict_int = 'MODULE_INTERACTION_GAIN'
    elif nsig['FLAT'] >= 1 or maxflat >= 4:
        verdict_int = 'NO_GAIN'
    else:
        verdict_int = 'NOT_ASSESSABLE'
    bb = list(csv.DictReader(open(C.ATT + '/att_bbj_tc_domains.csv')))
    bflag = np.array([int(r['bbj_tc_domain']) for r in bb], dtype=bool)
    sec = {}
    if bflag.sum() > 0:
        sub = dz[bflag]
        sec = dict(n_domains=int(bflag.sum()),
                   mean_delta=float(np.nanmean(sub)),
                   frac_delta_pos=float(np.nanmean(sub[np.isfinite(sub)] > 0)),
                   n_sig_ATT=int(np.nansum(q['ATT'][bflag] < FDR)),
                   n_sig_FLAT=int(np.nansum(q['FLAT'][bflag] < FDR)),
                   n_sig_FIXED=int(np.nansum(q['FIXED'][bflag] < FDR)),
                   n_sig_INT=int(np.nansum(q['INT'][bflag] < FDR)),
                   mean_delta_int=float(np.nanmean(dzi[bflag])),
                   max_abs_z_INT=float(np.nanmax(np.abs(z['INT'][bflag]))),
                   max_abs_z_FLAT=float(np.nanmax(np.abs(z['FLAT'][bflag]))),
                   max_abs_z_ATT=float(np.nanmax(np.abs(z['ATT'][bflag]))))
    rows = []
    for k in ('FLAT', 'FIXED', 'SIZE', 'ATT', 'INT', 'INT_product_only',
              'ATT_GENEHOLD', 'ATT_MAF001', 'FLAT_MAF001', 'FIXED_MAF001'):
        if k not in z:
            continue
        rows.append(dict(arm=k, n_sig_fdr05=nsig[k],
                         max_abs_z=float(np.nanmax(np.abs(z[k]))),
                         median_abs_z=float(np.nanmedian(np.abs(z[k]))),
                         mean_abs_z=float(np.nanmean(np.abs(z[k]))),
                         frac_abs_z_ge2=float(np.nanmean(np.abs(z[k]) >= 2)),
                         n_evaluable=int(np.isfinite(z[k]).sum())))
    with open(C.ATT + '/att_summary.csv', 'w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)
        f.write('\n')
        f.write('metric,value\n')
        for k, v in [('mean_delta_att_minus_flat', mean_delta),
                     ('median_delta', float(np.nanmedian(dz))),
                     ('frac_delta_pos', float(np.nanmean(dz[good] > 0))),
                     ('null_95pct_of_20_shuffle_means', n95),
                     ('null_95pct_pooled_20000_gene_deltas', p95),
                     ('null_max_of_20_shuffle_means',
                      float(null_means.max())),
                     ('sign_test_one_sided_p', sign_p),
                     ('n_delta_pos', npos), ('n_delta_nonzero', nnz),
                     ('cond1_mean_delta_gt_null95', int(c1)),
                     ('cond2_signtest_lt_0.025', int(c2)),
                     ('cond3_nATT_gt_nFLAT', int(c3)),
                     ('ref_MDE_z_ratio', RULES['ref_MDE']),
                     ('ref_B_crossfit_z_1', RULES['ref_B_crossfit_z'][0]),
                     ('ref_B_crossfit_z_2', RULES['ref_B_crossfit_z'][1]),
                     ('verdict', verdict),
                     ('mean_delta_int_minus_flat', mean_delta_int),
                     ('median_delta_int', float(np.nanmedian(dzi))),
                     ('frac_delta_int_pos', float(np.nanmean(dzi[gi] > 0))),
                     ('null_95pct_of_20_shuffle_means_int', n95i),
                     ('null_95pct_pooled_gene_deltas_int', p95i),
                     ('sign_test_int_one_sided_p', sign_pi),
                     ('cond1_int_mean_delta_gt_null95', int(i1)),
                     ('cond2_int_signtest_lt_0.025', int(i2)),
                     ('cond3_int_nINT_gt_nFLAT', int(i3)),
                     ('verdict_int', verdict_int)]:
            f.write('%s,%s\n' % (k, v))
    out = dict(amendment='A5', post_hoc=True, rules=RULES, verdict=verdict,
               conditions=dict(mean_delta=mean_delta, null_95=n95,
                               null_95_pooled=p95, sign_p=sign_p,
                               n_sig=nsig, cond1=bool(c1), cond2=bool(c2),
                               cond3=bool(c3), max_abs_z_FLAT=maxflat,
                               max_abs_z_any=maxall),
               verdict_int=verdict_int,
               conditions_int=dict(mean_delta=mean_delta_int, null_95=n95i,
                                   null_95_pooled=p95i, sign_p=sign_pi,
                                   n_sig_INT=nsig['INT'],
                                   n_sig_INT_product_only=nsig[
                                       'INT_product_only'],
                                   cond1=bool(i1), cond2=bool(i2),
                                   cond3=bool(i3)),
               secondary_bbj_tc=sec,
               judge_sha256=hashlib.sha256(
                   open(__file__, 'rb').read()).hexdigest())
    with open(C.ATT + '/att_judgement.json', 'w') as f:
        json.dump(out, f, indent=1, ensure_ascii=False, default=float)
    print(json.dumps(out, indent=1, ensure_ascii=False, default=float))

def _bh(p):
    p = np.asarray(p, dtype=np.float64)
    m = np.isfinite(p)
    out = np.full(len(p), np.nan)
    pp = p[m]
    n = len(pp)
    if n == 0:
        return out
    qq = pp * n / (np.argsort(np.argsort(pp)) + 1.0)
    run = np.inf
    o = np.empty(n)
    for i in np.argsort(-pp):
        run = min(run, qq[i])
        o[i] = run
    out[m] = np.minimum(o, 1.0)
    return out

if __name__ == '__main__':
    main()
