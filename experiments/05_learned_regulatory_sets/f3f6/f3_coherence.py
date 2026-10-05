import os
import sys
import json
import warnings
import numpy as np

import f3f6_common as C

warnings.filterwarnings('ignore')

BOOT = C.BOOT
AX_VAR = [('staar_cons', 'st_cons'), ('staar_epi_active', 'st_epi_active'),
          ('staar_epi_repr', 'st_epi_repr'),
          ('staar_epi_trans', 'st_epi_trans'), ('staar_tf', 'st_tf'),
          ('staar_linsight', 'st_linsight'), ('cadd_phred', 'cadd_phred'),
          ('gnomad_loeuf', 'gn_loeuf'), ('gnomad_pli', 'gn_pli'),
          ('gnomad_mis_z', 'gn_mis_z')]
AX_GRAM = ['re2g_profile_cos', 're2g_target_cos']
AX_HOMOG = ['staar_cage_prom_homog', 'staar_genehancer_homog']
AX_RE2GVAR = ['re2g_dist_tss']
AXES = ([a + '_varratio' for a, _ in AX_VAR] + AX_GRAM + AX_HOMOG
        + [a + '_varratio' for a in AX_RE2GVAR])
HIGHER_BETTER = dict([(a + '_varratio', False) for a, _ in AX_VAR]
                     + [(a, True) for a in AX_GRAM]
                     + [(a, True) for a in AX_HOMOG]
                     + [(a + '_varratio', False) for a in AX_RE2GVAR])
NA = len(AXES)

def _prep_domain(d, rows):
    v_list = []
    m_list = []
    names = []
    vcols = list(d['vset_cols'])
    for nm, col in AX_VAR:
        x = d['vset'][rows, vcols.index(col)].astype(np.float64)
        v_list.append(x)
        m_list.append(np.isfinite(x))
        names.append(nm + '_varratio')
    re2 = d['n_re2g'][rows] > 0
    dt = np.nanmean(np.where(np.isfinite(d['re2g_dist_tss'][rows]),
                             d['re2g_dist_tss'][rows], np.nan), axis=1)
    dt = np.where(re2, dt, np.nan)
    v_list.append(dt)
    m_list.append(np.isfinite(dt))
    names.append('re2g_dist_tss_varratio')
    hg = []
    for j, cn in enumerate(list(d['staar_txt_cols'])):
        b = d['staar_txt'][rows, j].astype(np.float64)
        seen = d['staar_txt_seen'][rows] > 0
        hg.append((b, seen))
    n = len(rows)
    S = d['re2g_score'][rows].astype(np.float64)
    ok = np.isfinite(S)
    Sf = np.where(ok, S, 0.0)
    cnt = ok.sum(1)
    mu = np.where(cnt > 0, Sf.sum(1) / np.maximum(cnt, 1), 0.0)
    U = np.where(ok, Sf - mu[:, None], 0.0)
    nrm = np.sqrt((U ** 2).sum(1))
    elig_p = re2 & (nrm > 0) & (cnt >= 2)
    U = np.where(elig_p[:, None], U / np.maximum(nrm, 1e-12)[:, None], 0.0)
    U = U.astype(np.float32)
    G_p = (U @ U.T).astype(np.float32)
    ptr = d['tgt_indptr']
    ind = d['tgt_indices']
    loc = {}
    per = []
    for i, r in enumerate(rows):
        t = ind[ptr[r]:ptr[r + 1]]
        per.append(t)
        for x in t:
            if x not in loc:
                loc[x] = len(loc)
    V = len(loc)
    if V > 0:
        T = np.zeros((n, V), dtype=np.float32)
        for i, t in enumerate(per):
            if len(t) == 0:
                continue
            T[i, [loc[x] for x in t]] = 1.0
        nn = np.sqrt((T ** 2).sum(1))
        elig_t = re2 & (nn > 0)
        T = np.where(elig_t[:, None], T / np.maximum(nn, 1e-12)[:, None], 0.0)
        G_t = (T.astype(np.float32) @ T.astype(np.float32).T).astype(np.float32)
    else:
        elig_t = np.zeros(n, dtype=bool)
        G_t = np.zeros((n, n), dtype=np.float32)
    Vv = np.array(v_list)
    Mm = np.array(m_list, dtype=np.float64)
    Vf = np.where(Mm > 0, np.nan_to_num(Vv, nan=0.0), 0.0)
    return dict(v=Vv, m=Mm, vf=Vf, vf2=Vf ** 2,
                names=names, homog=hg, n_vocab=V,
                G=[G_p, G_t], elig=[elig_p.astype(np.float64),
                                    elig_t.astype(np.float64)])

def _stats(P, Z):
    b = Z.shape[0]
    out = np.full((b, NA), np.nan)
    mean_out = np.full((b, NA), np.nan)
    nv = P['v'].shape[0]
    cnt = Z @ P['m'].T
    s = Z @ P['vf'].T
    s2 = Z @ P['vf2'].T
    with np.errstate(invalid='ignore', divide='ignore'):
        mu = s / cnt
        var = s2 / cnt - mu ** 2
    var = np.where(cnt >= 2, var, np.nan)
    mu = np.where(cnt >= 1, mu, np.nan)
    for j in range(nv):
        name = P['names'][j]
        k = AXES.index(name)
        out[:, k] = var[:, j]
        mean_out[:, k] = mu[:, j]
    Z32 = Z.astype(np.float32)
    for gi, nmg in enumerate(AX_GRAM):
        k = AXES.index(nmg)
        e = Z @ P['elig'][gi]
        q = ((Z32 @ P['G'][gi]) * Z32).sum(1).astype(np.float64)
        with np.errstate(invalid='ignore', divide='ignore'):
            val = (q - e) / (e * (e - 1.0))
        out[:, k] = np.where(e >= 2, val, np.nan)
        mean_out[:, k] = out[:, k]
    for hi, nmh in enumerate(AX_HOMOG):
        k = AXES.index(nmh)
        bvec, seen = P['homog'][hi]
        cs = Z @ seen.astype(np.float64)
        sb = Z @ (bvec * seen)
        with np.errstate(invalid='ignore', divide='ignore'):
            p = sb / cs
        out[:, k] = np.where(cs >= 2, 1.0 - 2.0 * p * (1.0 - p), np.nan)
        mean_out[:, k] = np.where(cs >= 1, p, np.nan)
    return out, mean_out

def run_scope(d, pool_of_domain, mod_lab, tag, null_kinds, out_prefix,
              domains=None):
    res = dict((nk, dict(num=np.zeros((BOOT, NA)), den=np.zeros((BOOT, NA)),
                         onum=np.zeros(NA), oden=np.zeros(NA),
                         mnum=np.zeros((BOOT, NA)), mden=np.zeros((BOOT, NA)),
                         omnum=np.zeros(NA), omden=np.zeros(NA),
                         nmod=0, ovl=[], bal=[]))
               for nk in null_kinds)
    detail = []
    ndeg = np.zeros(NA)
    for pid, rows in pool_of_domain:
        lab = mod_lab[rows]
        P = _prep_domain(d, rows)
        Zall = np.ones((1, len(rows)))
        pool_stat, pool_mean = _stats(P, Zall)
        pool_var = pool_stat[0]
        pool_sd = np.sqrt(np.where(pool_var > 0, pool_var, np.nan))
        ulab = [u for u in np.unique(lab) if u >= 0]
        strat = dict((nk, C.strata_for(d, rows, nk)) for nk in null_kinds)
        for u in ulab:
            mm = lab == u
            nm = int(mm.sum())
            if nm < 2:
                continue
            zobs = np.zeros((1, len(rows)))
            zobs[0, mm] = 1.0
            so, mo = _stats(P, zobs)
            so = so[0]
            mo = mo[0]
            for nk in null_kinds:
                idx = C.matched_null(strat[nk], mm, boot=BOOT,
                                     seed=C.det_seed(tag, nk, pid, int(u)))
                Z = np.zeros((BOOT, len(rows)))
                Z[np.arange(BOOT)[:, None], idx] = 1.0
                sn, mn = _stats(P, Z)
                R = res[nk]
                R['nmod'] += 1
                ovl = (Z[:, mm].sum(1) / float(nm)).mean()
                R['ovl'].append(ovl)
                R['bal'].append([
                    float(np.nanmean(d['maf'][rows][mm])
                          - np.nanmean((Z @ d['maf'][rows]) / nm)),
                    float(np.nanmean(d['r2'][rows][mm])
                          - np.nanmean((Z @ d['r2'][rows]) / nm)),
                    float(np.nanmean(d['dist_center'][rows][mm])
                          - np.nanmean((Z @ d['dist_center'][rows]) / nm)),
                    float((d['n_re2g'][rows][mm] > 0).mean()
                          - np.nanmean((Z @ (d['n_re2g'][rows] > 0)
                                        .astype(float)) / nm))])
                for k, ax in enumerate(AXES):
                    if HIGHER_BETTER[ax]:
                        vo, vn = so[k], sn[:, k]
                    else:
                        if not np.isfinite(pool_var[k]) or pool_var[k] <= 0:
                            ndeg[k] += 1
                            continue
                        vo = so[k] / pool_var[k]
                        vn = sn[:, k] / pool_var[k]
                    if np.isfinite(vo):
                        R['onum'][k] += vo * nm
                        R['oden'][k] += nm
                    okn = np.isfinite(vn)
                    R['num'][:, k] += np.where(okn, vn * nm, 0.0)
                    R['den'][:, k] += np.where(okn, nm, 0.0)
                    sd = np.nan if (ax in AX_GRAM or ax in AX_HOMOG) \
                        else pool_sd[k]
                    if np.isfinite(sd) and sd > 0:
                        eo = (mo[k] - pool_mean[0, k]) / sd
                        en = (mn[:, k] - pool_mean[0, k]) / sd
                        if np.isfinite(eo):
                            R['omnum'][k] += eo * nm
                            R['omden'][k] += nm
                        oke = np.isfinite(en)
                        R['mnum'][:, k] += np.where(oke, en * nm, 0.0)
                        R['mden'][:, k] += np.where(oke, nm, 0.0)
                if nk == null_kinds[0]:
                    detail.append([pid, int(u), nm, float(ovl)]
                                  + [float(so[k]) for k in range(NA)])
    return res, detail, ndeg

def _emit(res, tag, scope, ndeg, nmod_tot, rows_out, mrows_out):
    for nk, R in res.items():
        with np.errstate(invalid='ignore', divide='ignore'):
            oa = R['onum'] / R['oden']
            na = R['num'] / R['den']
            oma = R['omnum'] / R['omden']
            nma = R['mnum'] / R['mden']
        for k, ax in enumerate(AXES):
            hb = HIGHER_BETTER[ax]
            o = oa[k]
            nn = na[:, k]
            nn = nn[np.isfinite(nn)]
            if not np.isfinite(o) or len(nn) == 0:
                rows_out.append([tag, scope, 'null' + nk, ax,
                                 'higher' if hb else 'lower', '', '', '', '',
                                 '', '', 0, int(ndeg[k]), nmod_tot,
                                 float(np.mean(R['ovl'])) if R['ovl'] else ''])
                continue
            n_better = int((nn >= o).sum() if hb else (nn <= o).sum())
            p = (n_better + 1.0) / (len(nn) + 1.0)
            rows_out.append([tag, scope, 'null' + nk, ax,
                             'higher' if hb else 'lower', float(o),
                             float(np.mean(nn)), float(np.min(nn)),
                             float(np.max(nn)),
                             float(np.quantile(nn, 0.025)),
                             float(np.quantile(nn, 0.975)), n_better,
                             int(ndeg[k]), nmod_tot,
                             float(np.mean(R['ovl'])) if R['ovl'] else ''])
            om = oma[k]
            nm_ = nma[:, k]
            nm_ = nm_[np.isfinite(nm_)]
            if np.isfinite(om) and len(nm_) > 0:
                pe = (min((nm_ >= om).sum(), (nm_ <= om).sum()) + 1.0) \
                     / (len(nm_) + 1.0) * 2.0
                mrows_out.append([tag, scope, 'null' + nk, ax, float(om),
                                  float(np.mean(nm_)),
                                  float(np.quantile(nm_, 0.025)),
                                  float(np.quantile(nm_, 0.975)),
                                  min(1.0, float(pe)), nmod_tot])
        bal = np.array(R['bal']) if R['bal'] else np.zeros((1, 4))
        rows_out.append([tag, scope, 'null' + nk, '_BALANCE_maf_r2_dist_re2g',
                         '', float(np.nanmean(bal[:, 0])),
                         float(np.nanmean(bal[:, 1])),
                         float(np.nanmean(bal[:, 2])),
                         float(np.nanmean(bal[:, 3])), '', '', '', '',
                         nmod_tot, float(np.mean(R['ovl'])) if R['ovl'] else ''])

def main():
    os.makedirs(C.OUT, exist_ok=True)
    d = C.load()
    rows_out = []
    mrows_out = []
    detail_all = []

    lab = d['mod_spectral_eucl']
    di = d['dom_idx']
    nd = len(d['udom'])
    pools_all = [(int(g), np.where(di == g)[0]) for g in range(nd)]
    r701 = np.where(d['dom_forced'] == 0)[0]
    s701 = set(r701.tolist())
    for scope, pools in [('all1000', pools_all),
                         ('rand701', [p for p in pools_all if p[0] in s701])]:
        nks = ['A', 'B', 'C'] if scope == 'all1000' else ['B']
        res, det, ndeg = run_scope(d, pools, lab, 'modelC_' + scope, nks,
                                   'C')
        nmod = res[nks[0]]['nmod']
        _emit(res, 'modelC_spectral_eucl', scope, ndeg, nmod, rows_out,
              mrows_out)
        if scope == 'all1000':
            detail_all = det

    z = np.load(C.F1F2 + '/f2U_consensus_assign.npz', allow_pickle=True)
    uk = z['key']
    ud = z['domain']
    ul = z['spectral|eucl']
    rid = {}
    for i in range(len(d['key'])):
        rid[(d['key'][i], d['domain'][i])] = i
    rr = []
    ll = []
    for j in range(len(uk)):
        i = rid.get((uk[j], ud[j]))
        if i is not None:
            rr.append(i)
            ll.append(int(ul[j]))
    rr = np.array(rr)
    labU = np.full(len(d['key']), -1, dtype=np.int32)
    labU[rr] = np.array(ll, dtype=np.int32)
    resU, detU, ndegU = run_scope(d, [(-1, rr)], labU, 'modelU', ['B'], 'U')
    _emit(resU, 'modelU_spectral_eucl', 'subsample8000', ndegU,
          resU['B']['nmod'], rows_out, mrows_out)

    hdr = ['model', 'scope', 'null', 'axis', 'coherent_direction', 'observed',
           'null_mean', 'null_min', 'null_max', 'null_q025', 'null_q975',
           'n_null_at_least_as_coherent', 'n_modules_degenerate_axis',
           'n_modules_evaluated', 'mean_null_module_overlap']
    with open(C.OUT + '/biological_coherence.csv', 'w') as fh:
        fh.write(','.join(hdr) + '\n')
        for r in rows_out:
            fh.write(','.join([str(x) for x in r]) + '\n')
    hdr2 = ['model', 'scope', 'null', 'axis', 'observed_std_mean_diff',
            'null_mean', 'null_q025', 'null_q975', 'p_two_sided',
            'n_modules_evaluated']
    with open(C.OUT + '/module_enrichment.csv', 'w') as fh:
        fh.write(','.join(hdr2) + '\n')
        for r in mrows_out:
            fh.write(','.join([str(x) for x in r]) + '\n')
    import gzip as _gz
    with _gz.open(C.OUT + '/f3_module_detail.csv.gz', 'wt') as fh:
        fh.write(','.join(['domain_idx', 'module', 'n', 'null_overlap']
                          + AXES) + '\n')
        for r in detail_all:
            fh.write(','.join([str(x) for x in r]) + '\n')
    with open(C.OUT + '/f3.done', 'w') as fh:
        json.dump(dict(axes=AXES, boot=BOOT,
                       n_rows_C=int(len(d['key'])), n_rows_U=int(len(rr))),
                  fh, indent=1)

if __name__ == '__main__':
    main()
