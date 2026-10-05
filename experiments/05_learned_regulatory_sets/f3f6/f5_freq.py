import os
import json
import numpy as np
from scipy.stats import kruskal

import f3f6_common as C
import f3_coherence as F3

BOOT = C.BOOT
RULE = ('F5 PASS <=> observed aggregate MI(module, MAF decile) <= 95th '
        'percentile of 1,000 within-domain label permutations')

def mi_from_counts(T):
    N = T.sum()
    if N <= 0:
        return np.nan, np.nan
    P = T / N
    pr = P.sum(1, keepdims=True)
    pc = P.sum(0, keepdims=True)
    with np.errstate(divide='ignore', invalid='ignore'):
        L = np.where(P > 0, P * np.log2(P / (pr * pc)), 0.0)
    mi = float(L.sum())
    hr = float(-(pr[pr > 0] * np.log2(pr[pr > 0])).sum())
    hc = float(-(pc[pc > 0] * np.log2(pc[pc > 0])).sum())
    den = 0.5 * (hr + hc)
    return mi, (mi / den if den > 0 else np.nan)

def main():
    d = C.load()
    lab = d['mod_spectral_eucl']
    di = d['dom_idx']
    nd = len(d['udom'])
    r701 = set(np.where(d['dom_forced'] == 0)[0].tolist())

    out = []
    per_dom = []
    for scope in ['all1000', 'rand701']:
        doms = [g for g in range(nd) if (scope == 'all1000' or g in r701)]
        num_o = 0.0
        num_n = np.zeros(BOOT)
        den = 0.0
        numN_o = 0.0
        numN_n = np.zeros(BOOT)
        denN = 0.0
        kw_p = []
        for g in doms:
            rid = np.where(di == g)[0]
            lk = lab[rid]
            u = np.array([x for x in np.unique(lk) if x >= 0])
            if len(u) < 2:
                continue
            v = d['maf'][rid]
            ok = np.isfinite(v)
            if ok.sum() < 10:
                continue
            q = np.quantile(v[ok], np.arange(1, 5) / 5.0)
            b = np.digitize(v, q)
            m = np.isin(lk, u) & ok
            li = np.searchsorted(u, lk[m])
            bi = b[m]
            nb = int(bi.max()) + 1
            T = np.zeros((len(u), nb))
            np.add.at(T, (li, bi), 1.0)
            mi, nmi = mi_from_counts(T)
            w = float(m.sum())
            num_o += mi * w
            numN_o += nmi * w
            den += w
            rng = np.random.default_rng(C.det_seed('f5', int(g)))
            mis = np.empty(BOOT)
            nmis = np.empty(BOOT)
            for t in range(BOOT):
                lp = li[rng.permutation(len(li))]
                Tp = np.zeros((len(u), nb))
                np.add.at(Tp, (lp, bi), 1.0)
                a, bb = mi_from_counts(Tp)
                mis[t] = a
                nmis[t] = bb
            num_n += mis * w
            numN_n += nmis * w
            denN += w
            try:
                gg = [v[m][li == j] for j in range(len(u))]
                gg = [x for x in gg if len(x) > 0]
                kp = float(kruskal(*gg)[1]) if len(gg) >= 2 else np.nan
            except Exception:
                kp = np.nan
            kw_p.append(kp)
            if scope == 'all1000':
                per_dom.append([int(g), int(d['dom_forced'][g]), len(u),
                                int(w), mi, nmi, float(np.mean(mis)),
                                float(np.quantile(mis, 0.95)), kp])
        mo = num_o / den
        mn = num_n / den
        q95 = float(np.quantile(mn, 0.95))
        nmo = numN_o / denN
        nmn = numN_n / denN
        kw = np.array(kw_p, dtype=float)
        kwok = np.isfinite(kw)
        out.append(['modelC_spectral_eucl', scope, 'MI_bits_weighted',
                    float(mo), float(np.mean(mn)), q95,
                    float(np.quantile(mn, 0.025)),
                    float(np.quantile(mn, 0.975)),
                    int((mn >= mo).sum()), 'MI<=null_q95',
                    'PASS' if mo <= q95 else 'FAIL', len(doms)])
        out.append(['modelC_spectral_eucl', scope, 'NMI_weighted', float(nmo),
                    float(np.mean(nmn)), float(np.quantile(nmn, 0.95)),
                    float(np.quantile(nmn, 0.025)),
                    float(np.quantile(nmn, 0.975)),
                    int((nmn >= nmo).sum()), 'reference', '', len(doms)])
        out.append(['modelC_spectral_eucl', scope, 'KW_frac_p_lt_0.05',
                    float((kw[kwok] < 0.05).mean()) if kwok.any() else '',
                    0.05, '', '', '', '', 'reference', '',
                    int(kwok.sum())])

    pools = [(int(g), np.where(di == g)[0]) for g in range(nd)]
    resM, _, ndegM = F3.run_scope(d, pools, lab, 'f5nullM', ['M'], 'M')
    rows_c = []
    mrows_c = []
    F3._emit(resM, 'modelC_spectral_eucl', 'all1000_nullM', ndegM,
             resM['M']['nmod'], rows_c, mrows_c)

    hdr = ['model', 'scope', 'statistic', 'observed', 'null_mean',
           'null_q95', 'null_q025', 'null_q975', 'n_null_ge_obs', 'rule',
           'verdict', 'n_domains']
    with open(C.OUT + '/frequency_independence.csv', 'w') as fh:
        fh.write(','.join(hdr) + '\n')
        for r in out:
            fh.write(','.join([str(x) for x in r]) + '\n')
        fh.write('\n')
        fh.write('# Null M (has_re2g x MAF decile exact match) V-set coherence'
                 '\n')
        fh.write(','.join(['model', 'scope', 'null', 'axis',
                           'coherent_direction', 'observed', 'null_mean',
                           'null_min', 'null_max', 'null_q025', 'null_q975',
                           'n_null_at_least_as_coherent',
                           'n_modules_degenerate_axis',
                           'n_modules_evaluated',
                           'mean_null_module_overlap']) + '\n')
        for r in rows_c:
            fh.write(','.join([str(x) for x in r]) + '\n')
    with open(C.OUT + '/f5_per_domain.csv', 'w') as fh:
        fh.write(','.join(['domain_idx', 'forced', 'k_modules', 'n_rows',
                           'mi_bits', 'nmi', 'null_mi_mean', 'null_mi_q95',
                           'kruskal_p']) + '\n')
        for r in per_dom:
            fh.write(','.join([str(x) for x in r]) + '\n')
    v = [r for r in out if r[2] == 'MI_bits_weighted' and r[1] == 'all1000']
    with open(C.OUT + '/f5_judgement.json', 'w') as fh:
        json.dump(dict(rule=RULE, verdict=(v[0][10] if v else 'NA'),
                       observed=(v[0][3] if v else None),
                       null_q95=(v[0][5] if v else None), boot=BOOT), fh,
                  indent=1)
    with open(C.OUT + '/f5.done', 'w') as fh:
        fh.write((v[0][10] if v else 'NA') + '\n')

if __name__ == '__main__':
    main()
