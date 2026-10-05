#!/usr/bin/env python3
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import csv, json, os, math, statistics as stt
import numpy as np
W = _config_path('${PROJECT_ROOT}/work/run_trackA')
def f_(x): return float('nan') if x in (None, 'NA', '') else float(x)
def read(p): return list(csv.DictReader(open(p), delimiter='\t')) if os.path.exists(p) else []
def sizes(s): return sorted(int(x) for x in s.split(';') if x) if s else []
cov = {r['region']: r for r in read(f'{W}/phi/phi_coverage_all.tsv')}; diag = {r['region']: r for r in read(f'{W}/phi/phi_diag.tsv')}
regs = [l.split('\t') for l in open(f'{W}/regions_by_size_v2.tsv').read().strip().split('\n')]
def pcts(phi_cs, phi_mp, phi_en, pc, pm, pe):
    n = len(pm)
    return dict(cs=0.0 if math.isnan(phi_cs) else sum(1 for x in pc if (not math.isnan(x)) and x > phi_cs) / n, mp=sum(1 for x in pm if x < phi_mp) / n, en=sum(1 for x in pe if x > phi_en) / n)
def rank_corr(a, b):
    from scipy.stats import spearmanr; return float(spearmanr(a, b).correlation)
out = []
for t, rid, nfm in regs:
    os_ = read(_config_path(f'${{PROJECT_ROOT}}/work/run_ourfm/fm/{t}/{rid}.summary.tsv')); os_ = os_[0] if os_ else {}
    u_sizes = sizes(os_.get('cs_sizes', ''))
    row = dict(trait=t, region=rid, n_fm=int(nfm), u_converged=os_.get('converged'), u_n_cs=len(u_sizes), u_mean_cs=(sum(u_sizes) / len(u_sizes)) if u_sizes else float('nan'), u_max_pip=f_(os_.get('max_pip')))
    c = cov.get(rid); d = diag.get(rid)
    if c: row.update(frac_median=float(c['frac_median']), phi_ratio_max_min=float(c['phi_ratio_max_min']), phi_cv=float(c['phi_cv']))
    if d: row.update(rho_phi_dist=f_(d['rho_phi_dist']), rho_phi_logdist=f_(d['rho_phi_logdist']), frac_dist_na=f_(d['frac_dist_na']), phi_p99_p01=f_(d['phi_p99_p01']), phi_top1pct_mass=f_(d['phi_top1pct_mass']))
    g = {r['arm']: r for r in read(f'{W}/fm/{rid}.gate.arms.tsv')}.get('gate_uniform') or {r['arm']: r for r in read(f'{W}/fm/{rid}.arms.tsv')}.get('gate_uniform')
    rwa = {r['arm']: r for r in read(f'{W}/rw/{rid}.rw.arms.tsv')}
    if g: row.update(gate_max_abs_dpip=f_(g['gate_max_abs_dpip']), gate_cs_equal=str(sizes(g['cs_sizes']) == u_sizes).upper())
    elif 'uniform_refit' in rwa: row.update(gate_max_abs_dpip=f_(rwa['uniform_refit']['gate_max_abs_dpip']), gate_cs_equal=str(sizes(rwa['uniform_refit']['cs_sizes']) == u_sizes).upper(), gate_src='rw')
    if 'phi_rwA' in rwa:
        ur = rwa['uniform_refit']; row.update(u_n_cs_pure=int(ur['n_cs']), u_n_cs_all=int(ur['n_cs_all']), u_entropy=f_(ur['entropy']), u_cs_sizes=ur['cs_sizes'], u_purity=ur['cs_min_abs_corr'])
        perms = read(f'{W}/rw/{rid}.rw.perms.tsv'); row['n_perm_rw'] = len(perms)
        for V in ('A', 'B'):
            p = rwa[f'phi_rw{V}']; row[f'{V}_n_cs'] = int(p['n_cs']); row[f'{V}_mean_cs'] = f_(p['mean_cs_size']); row[f'{V}_max_pip'] = f_(p['max_pip']); row[f'{V}_entropy'] = f_(p['entropy']); row[f'{V}_n_eff_rw'] = p['n_eff_reweighted']; row[f'{V}_cs_sizes'] = p['cs_sizes']
            pc = [f_(r[f'{V}_mean_cs']) for r in perms]; pm = [f_(r[f'{V}_max_pip']) for r in perms]; pe = [f_(r[f'{V}_entropy']) for r in perms]
            q = pcts(row[f'{V}_mean_cs'], row[f'{V}_max_pip'], row[f'{V}_entropy'], pc, pm, pe)
            row[f'{V}_pct_cs_shrink'] = q['cs']; row[f'{V}_pct_maxpip_up'] = q['mp']; row[f'{V}_pct_entropy_down'] = q['en']; row[f'{V}_exceeds'] = int(q['cs'] >= 0.95 or q['mp'] >= 0.95)
            pcv = [x for x in pc if not math.isnan(x)]; row[f'{V}_perm_mean_cs_med'] = stt.median(pcv) if pcv else float('nan'); row[f'{V}_perm_max_pip_med'] = stt.median(pm); row[f'{V}_perm_max_pip_p95'] = float(np.percentile(pm, 95))
        row['rw_status'] = 'done'
    else: row['rw_status'] = 'pending'
    a = read(f'{W}/fm/{rid}.arms.tsv'); da = {r['arm']: r for r in a}; perms20 = [r for r in a if r['arm'].startswith('perm')]
    if 'phi' in da and perms20:
        p = da['phi']; row.update(R_n_cs=int(p['n_cs']), R_mean_cs=f_(p['mean_cs_size']), R_max_pip=f_(p['max_pip']), R_entropy=f_(p['entropy']), R_cs_sizes=p['cs_sizes'], R_converged=p['converged'], n_perm_refit=len(perms20))
        q = pcts(row['R_mean_cs'], row['R_max_pip'], row['R_entropy'], [f_(r['mean_cs_size']) for r in perms20], [f_(r['max_pip']) for r in perms20], [f_(r['entropy']) for r in perms20])
        row['R_pct_cs_shrink'] = q['cs']; row['R_pct_maxpip_up'] = q['mp']; row['R_pct_entropy_down'] = q['en']; row['R_exceeds'] = int(q['cs'] >= 0.95 or q['mp'] >= 0.95); row['refit_status'] = 'done'
        if os.path.exists(f'{W}/fm/{rid}.phi.pip.tsv') and os.path.exists(f'{W}/rw/{rid}.rwA.pip.tsv'):
            pr = read(f'{W}/fm/{rid}.phi.pip.tsv'); pa = read(f'{W}/rw/{rid}.rwA.pip.tsv'); pb = read(f'{W}/rw/{rid}.rwB.pip.tsv')
            assert [r['variant_hg19'] for r in pr] == [r['variant_hg19'] for r in pa]
            x = np.array([float(r['pip']) for r in pr]); ya = np.array([float(r['pip']) for r in pa]); yb = np.array([float(r['pip']) for r in pb])
            row['corr_pip_R_vs_A'] = float(np.corrcoef(x, ya)[0, 1]); row['spearman_pip_R_vs_A'] = rank_corr(x, ya); row['corr_pip_R_vs_B'] = float(np.corrcoef(x, yb)[0, 1]); row['max_abs_dpip_R_vs_A'] = float(np.abs(x - ya).max())
            row['cs_equal_R_vs_A'] = str(sizes(row['R_cs_sizes']) == sizes(row['A_cs_sizes'])).upper() if 'A_cs_sizes' in row else 'NA'
            if 'A_exceeds' in row: row['verdict_agree_R_vs_A'] = str(row['R_exceeds'] == row['A_exceeds']).upper()
    else: row['refit_status'] = 'pending'
    out.append(row)
cols = []
for r in out:
    for k in r:
        if k not in cols: cols.append(k)
with open(f'{W}/trackA_region_metrics_v2.tsv', 'w', newline='') as o:
    w = csv.DictWriter(o, fieldnames=cols, delimiter='\t', restval='NA'); w.writeheader()
    for r in out: w.writerow({k: (f'{v:.5g}' if isinstance(v, float) else v) for k, v in r.items()})
def verdict(rows, key):
    rows = [r for r in rows if f'{key}_exceeds' in r]
    if not rows: return dict(n=0)
    ne = sum(r[f'{key}_exceeds'] for r in rows); n = len(rows)
    return dict(n=n, n_exceeds=ne, frac=ne / n, n_cs_shrink=sum(1 for r in rows if r[f'{key}_pct_cs_shrink'] >= .95), n_maxpip=sum(1 for r in rows if r[f'{key}_pct_maxpip_up'] >= .95), n_entropy=sum(1 for r in rows if r[f'{key}_pct_entropy_down'] >= .95),
                med_pct_cs=stt.median(r[f'{key}_pct_cs_shrink'] for r in rows), med_pct_mp=stt.median(r[f'{key}_pct_maxpip_up'] for r in rows), med_pct_en=stt.median(r[f'{key}_pct_entropy_down'] for r in rows),
                verdict=('유효' if ne / n > .5 else '무효(기준 미달)') + ('' if n == 81 else ' — 잠정'), by_trait={t: [sum(1 for r in rows if r['trait'] == t), sum(r[f'{key}_exceeds'] for r in rows if r['trait'] == t)] for t in ['tchl', 'htn', 'dm', 'lip']})
rwd = [r for r in out if r['rw_status'] == 'done']
summ = dict(n_regions=len(out), rw_done=len(rwd), refit_done=sum(1 for r in out if r['refit_status'] == 'done'),
            gate=dict(n=sum(1 for r in out if 'gate_max_abs_dpip' in r), n_pass=sum(1 for r in out if 'gate_max_abs_dpip' in r and r['gate_max_abs_dpip'] <= 1e-6 and r['gate_cs_equal'] == 'TRUE'), max_dpip=max([r['gate_max_abs_dpip'] for r in out if 'gate_max_abs_dpip' in r], default=None)),
            A_ems_pure=verdict(out, 'A'), B_all_effects=verdict(out, 'B'), R_refit_perm20=verdict(out, 'R'),
            pure_cs=dict(regions_with_pure_cs=sum(1 for r in rwd if r['u_n_cs_pure'] > 0), regions_zero_pure_cs=sum(1 for r in rwd if r['u_n_cs_pure'] == 0), total_pure_cs=sum(r['u_n_cs_pure'] for r in rwd), total_cs_all=sum(r['u_n_cs_all'] for r in rwd),
                         A_verdict_among_pure_regions=verdict([r for r in rwd if r['u_n_cs_pure'] > 0], 'A'), A_verdict_among_zero_pure=verdict([r for r in rwd if r['u_n_cs_pure'] == 0], 'A'),
                         R_verdict_among_pure_regions=verdict([r for r in out if r.get('u_n_cs_pure', r['u_n_cs']) > 0], 'R'), R_verdict_among_zero_pure=verdict([r for r in out if r.get('u_n_cs_pure', r['u_n_cs']) == 0], 'R')),
            R_vs_A=dict(n=sum(1 for r in out if 'corr_pip_R_vs_A' in r), corr_pip_median=stt.median([r['corr_pip_R_vs_A'] for r in out if 'corr_pip_R_vs_A' in r]) if any('corr_pip_R_vs_A' in r for r in out) else None,
                        corr_pip_min=min([r['corr_pip_R_vs_A'] for r in out if 'corr_pip_R_vs_A' in r], default=None), cs_equal=sum(1 for r in out if r.get('cs_equal_R_vs_A') == 'TRUE'), verdict_agree=sum(1 for r in out if r.get('verdict_agree_R_vs_A') == 'TRUE'), verdict_compared=sum(1 for r in out if 'verdict_agree_R_vs_A' in r)),
            dist_tss=dict(n=sum(1 for r in out if 'rho_phi_dist' in r and not math.isnan(r['rho_phi_dist'])), rho_median=stt.median([r['rho_phi_dist'] for r in out if 'rho_phi_dist' in r and not math.isnan(r['rho_phi_dist'])]) if any('rho_phi_dist' in r for r in out) else None,
                          rho_min=min([r['rho_phi_dist'] for r in out if 'rho_phi_dist' in r and not math.isnan(r['rho_phi_dist'])], default=None), rho_max=max([r['rho_phi_dist'] for r in out if 'rho_phi_dist' in r and not math.isnan(r['rho_phi_dist'])], default=None)))
def chars(rows):
    if not rows: return {}
    def med(k): v = [r[k] for r in rows if k in r and not (isinstance(r[k], float) and math.isnan(r[k]))]; return stt.median(v) if v else None
    return dict(n=len(rows), n_fm_med=med('n_fm'), u_n_cs_pure_med=med('u_n_cs_pure'), u_n_cs_all_med=med('u_n_cs_all'), frac_pure=(sum(r.get('u_n_cs_pure', 0) for r in rows) / max(1, sum(r.get('u_n_cs_all', 0) for r in rows))), u_mean_cs_med=med('u_mean_cs'), u_max_pip_med=med('u_max_pip'),
                phi_ratio_max_min_med=med('phi_ratio_max_min'), phi_cv_med=med('phi_cv'), phi_p99_p01_med=med('phi_p99_p01'), rho_phi_dist_med=med('rho_phi_dist'), frac_median_med=med('frac_median'), u_converged_true=sum(1 for r in rows if r['u_converged'] == 'TRUE'))
summ['working_vs_not_A'] = dict(working=chars([r for r in rwd if r.get('A_exceeds') == 1]), not_working=chars([r for r in rwd if r.get('A_exceeds') == 0]))
summ['working_vs_not_R'] = dict(working=chars([r for r in out if r.get('R_exceeds') == 1]), not_working=chars([r for r in out if r.get('R_exceeds') == 0]))
json.dump(summ, open(f'{W}/trackA_summary_v2.json', 'w'), indent=1, ensure_ascii=False)
print(json.dumps({k: summ[k] for k in ['rw_done', 'refit_done', 'gate', 'A_ems_pure', 'B_all_effects', 'R_refit_perm20', 'R_vs_A', 'dist_tss']}, ensure_ascii=False))
