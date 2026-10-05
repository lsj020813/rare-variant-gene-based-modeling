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


import glob, csv, json, os, math, statistics as stt
W = _config_path('${PROJECT_ROOT}/work/run_trackA')
def f_(x):
    return float('nan') if x in (None, 'NA', '') else float(x)
def read(p): return list(csv.DictReader(open(p), delimiter='\t')) if os.path.exists(p) else []
cov = {r['region']: r for r in read(f'{W}/phi/phi_coverage_all.tsv')}
def orig_summary(t, rid):
    p = _config_path(f'${{PROJECT_ROOT}}/work/run_ourfm/fm/{t}/{rid}.summary.tsv'); r = read(p)
    return r[0] if r else {}
def sizes(s): return sorted(int(x) for x in s.split(';') if x) if s else []
regs = [l.split('\t') for l in open(f'{W}/regions_by_size.tsv').read().strip().split('\n')]
out = []; gate_rows = []
for t, rid, nfm in regs:
    g = {r['arm']: r for r in read(f'{W}/fm/{rid}.gate.arms.tsv')}
    a = read(f'{W}/fm/{rid}.arms.tsv'); d = {r['arm']: r for r in a}; perms = [r for r in a if r['arm'].startswith('perm')]
    row = dict(trait=t, region=rid, n_fm=int(nfm), status='pending')
    os_ = orig_summary(t, rid); u_sizes = sizes(os_.get('cs_sizes', ''))
    row.update(u_n_cs=len(u_sizes), u_mean_cs=(sum(u_sizes) / len(u_sizes)) if u_sizes else float('nan'), u_converged=os_.get('converged'), u_niter=os_.get('niter'))
    if g:
        gu = g.get('gate_uniform', {})
        row.update(gate_max_abs_dpip=f_(gu.get('gate_max_abs_dpip')), gate_cs_equal=str(sizes(gu.get('cs_sizes', '')) == u_sizes).upper(), gate_cs_equal_pipfile=gu.get('gate_cs_equal'), gate_niter=gu.get('niter'), gate_converged=gu.get('converged'))
        gate_rows.append(row)
    c = cov.get(rid)
    if c: row.update(n_bbj=int(c['n_bbj']), n_trA=int(c['n_trA']), n_median=int(c['n_median']), frac_median=float(c['frac_median']), phi_ratio_max_min=float(c['phi_ratio_max_min']), phi_cv=float(c['phi_cv']))
    if 'phi' in d and len(perms) >= 1:
        u = d.get('uniform_ref') or g.get('uniform_ref'); p = d['phi']
        row['u_max_pip'] = f_(u['max_pip']); row['u_entropy'] = f_(u['entropy'])
        row['phi_n_cs'] = int(p['n_cs']); row['phi_mean_cs'] = f_(p['mean_cs_size']); row['phi_max_pip'] = f_(p['max_pip']); row['phi_entropy'] = f_(p['entropy']); row['phi_cs_sizes'] = p['cs_sizes']; row['u_cs_sizes'] = os_.get('cs_sizes', '')
        row['phi_converged'] = p['converged']; row['phi_niter'] = p['niter']; row['n_perm'] = len(perms)
        pc = [f_(r['mean_cs_size']) for r in perms]; pm = [f_(r['max_pip']) for r in perms]; pe = [f_(r['entropy']) for r in perms]; pn = [int(r['n_cs']) for r in perms]
        row['perm_mean_cs_med'] = stt.median([x for x in pc if not math.isnan(x)]) if any(not math.isnan(x) for x in pc) else float('nan')
        row['perm_mean_cs_p05'] = sorted([x for x in pc if not math.isnan(x)])[max(0, int(0.05 * sum(1 for x in pc if not math.isnan(x))) - 1)] if any(not math.isnan(x) for x in pc) else float('nan')
        row['perm_max_pip_med'] = stt.median(pm); row['perm_max_pip_p95'] = sorted(pm)[min(len(pm) - 1, int(math.ceil(0.95 * len(pm))) - 1)]
        row['perm_entropy_med'] = stt.median(pe); row['perm_n_cs_med'] = stt.median(pn)
        n = len(perms)
        row['pct_cs_shrink'] = 0.0 if math.isnan(row['phi_mean_cs']) else sum(1 for x in pc if (not math.isnan(x)) and x > row['phi_mean_cs']) / n
        row['pct_maxpip_up'] = sum(1 for x in pm if x < row['phi_max_pip']) / n
        row['pct_entropy_down'] = sum(1 for x in pe if x > row['phi_entropy']) / n
        row['pct_ncs_up'] = sum(1 for x in pn if x < row['phi_n_cs']) / n
        row['d_mean_cs_vs_u'] = row['phi_mean_cs'] - row['u_mean_cs']; row['d_max_pip_vs_u'] = row['phi_max_pip'] - row['u_max_pip']; row['d_entropy_vs_u'] = row['phi_entropy'] - row['u_entropy']; row['d_n_cs_vs_u'] = row['phi_n_cs'] - row['u_n_cs']
        row['exceeds_p95'] = int(row['pct_cs_shrink'] >= 0.95 or row['pct_maxpip_up'] >= 0.95)
        row['status'] = 'done' if n >= 20 else f'partial_perm{n}'
        row['elapsed_total_s'] = round(sum(f_(r['elapsed_s']) for r in a if not math.isnan(f_(r['elapsed_s']))))
    elif os.path.exists(f'{W}/fm/{rid}.skip'): row['status'] = 'skip'
    out.append(row)
cols = []
for r in out:
    for k in r:
        if k not in cols: cols.append(k)
with open(f'{W}/trackA_region_metrics.tsv', 'w', newline='') as o:
    w = csv.DictWriter(o, fieldnames=cols, delimiter='\t', restval='NA'); w.writeheader()
    for r in out: w.writerow({k: (f'{v:.5g}' if isinstance(v, float) else v) for k, v in r.items()})
done = [r for r in out if r['status'] == 'done']
gd = [r for r in gate_rows if not math.isnan(r.get('gate_max_abs_dpip', float('nan')))]
summ = dict(n_regions=len(out), n_gate_done=len(gd), gate_pass=sum(1 for r in gd if r['gate_max_abs_dpip'] <= 1e-6 and r['gate_cs_equal'] == 'TRUE'),
            gate_max_abs_dpip_max=max([r['gate_max_abs_dpip'] for r in gd], default=float('nan')),
            n_done=len(done), n_partial=sum(1 for r in out if r['status'].startswith('partial')), n_skip=sum(1 for r in out if r['status'] == 'skip'),
            n_exceeds_p95=sum(r['exceeds_p95'] for r in done), frac_exceeds=(sum(r['exceeds_p95'] for r in done) / len(done)) if done else None,
            n_exceeds_cs=sum(1 for r in done if r['pct_cs_shrink'] >= 0.95), n_exceeds_maxpip=sum(1 for r in done if r['pct_maxpip_up'] >= 0.95),
            n_exceeds_entropy=sum(1 for r in done if r['pct_entropy_down'] >= 0.95),
            median_pct_cs_shrink=stt.median([r['pct_cs_shrink'] for r in done]) if done else None, median_pct_maxpip_up=stt.median([r['pct_maxpip_up'] for r in done]) if done else None,
            median_pct_entropy_down=stt.median([r['pct_entropy_down'] for r in done]) if done else None,
            n_phi_cs_fewer=sum(1 for r in done if r['d_n_cs_vs_u'] < 0), n_phi_cs_more=sum(1 for r in done if r['d_n_cs_vs_u'] > 0),
            n_phi_meancs_smaller=sum(1 for r in done if not math.isnan(r['d_mean_cs_vs_u']) and r['d_mean_cs_vs_u'] < 0), n_phi_maxpip_higher=sum(1 for r in done if r['d_max_pip_vs_u'] > 0),
            n_phi_entropy_lower=sum(1 for r in done if r['d_entropy_vs_u'] < 0),
            verdict=('유효' if done and sum(r['exceeds_p95'] for r in done) / len(done) > 0.5 else '무효(기준 미달)') + ('' if len(done) == len(out) else ' — 잠정(부분 완료)'),
            by_trait={t: dict(n_done=sum(1 for r in done if r['trait'] == t), n_exceeds=sum(r['exceeds_p95'] for r in done if r['trait'] == t)) for t in ['tchl', 'htn', 'dm', 'lip']})
json.dump(summ, open(f'{W}/trackA_summary.json', 'w'), indent=1, ensure_ascii=False)
print(json.dumps(summ, ensure_ascii=False))
