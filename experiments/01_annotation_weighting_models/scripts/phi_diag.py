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


import sys, os, glob, gzip, math
sys.path.insert(0, _config_path('${PROJECT_ROOT}/work/run_band15/model_v10_out')); import l1_train_v10 as L; L.load_libraries(4)
import numpy as np; from scipy.stats import spearmanr
W = _config_path('${PROJECT_ROOT}/work/run_ourfm'); OUT = _config_path('${PROJECT_ROOT}/work/run_trackA')
def num(v): return float('nan') if v in ('', 'NA') else float(v)
def swap(k): a = k.split(':'); return f'{a[0]}:{a[1]}:{a[3]}:{a[2]}'
trA = {}
with gzip.open(f'{OUT}/annot_trA/trA_annot_missing.tsv.gz', 'rt') as fh:
    hdr = fh.readline().rstrip('\n').split('\t'); ki = hdr.index('variant_hg19'); di = hdr.index('dist_tss'); ti = hdr.index('in_tss3kb'); na = hdr.index('t1_na')
    for line in fh: a = line.rstrip('\n').split('\t'); trA[a[ki]] = (num(a[di]), num(a[ti]), num(a[na]))
model_cols = __import__('json').load(open(_config_path('${PROJECT_ROOT}/work/run_band15/model_v10_out/e6_runs/primary_f0/phi.json')))['design']['selected_columns']
rows = []
for t in ['dm', 'htn', 'lip', 'tchl']:
    for f in sorted(glob.glob(f'{W}/fm/{t}/*.ld.vars')):
        rid = os.path.basename(f)[:-len('.ld.vars')]
        if '.mi1000' in rid: continue
        bbj = {}
        bf = f'{OUT}/annot/{rid}.bbj.tsv'
        if os.path.exists(bf):
          with open(bf) as fh:
            hdr = fh.readline().rstrip('\n').split('\t'); ki = hdr.index('variant_hg19'); di = hdr.index('dist_tss'); ti = hdr.index('in_tss3kb'); na = hdr.index('t1_na')
            for line in fh: a = line.rstrip('\n').split('\t'); k = a[ki] if a[0] == '0' else swap(a[ki]); bbj[k] = (num(a[di]), num(a[ti]), num(a[na]))
        ph = {}; 
        with open(f'{OUT}/phi/{rid}.phi.tsv') as fh:
            next(fh)
            for line in fh: k, p, s_ = line.rstrip('\n').split('\t'); ph[k] = float(p)
        phi = []; dist = []; tss3 = []; t1 = []
        for k, p in ph.items():
            r = bbj.get(k) or trA.get(k)
            if r is None: continue
            phi.append(p); dist.append(r[0]); tss3.append(r[1]); t1.append(r[2])
        phi = np.array(phi); dist = np.array(dist); tss3 = np.array(tss3); t1 = np.array(t1)
        ok = np.isfinite(dist)
        rho = spearmanr(phi[ok], dist[ok]).correlation if ok.sum() > 10 else float('nan')
        rho_log = spearmanr(phi[ok], np.log1p(dist[ok])).correlation if ok.sum() > 10 else float('nan')
        d2 = np.where(ok, dist, (np.nanmax(dist) if ok.any() else 0) + 1); rho_all = spearmanr(phi, d2).correlation
        rows.append(dict(trait=t, region=rid, n=len(phi), frac_dist_na=round(float((~ok).mean()), 4), rho_phi_dist=round(float(rho), 4), rho_phi_logdist=round(float(rho_log), 4), rho_phi_dist_naAsFar=round(float(rho_all), 4),
                         phi_mean_tss3kb=round(float(phi[tss3 == 1].mean()), 5) if (tss3 == 1).any() else float('nan'), phi_mean_not_tss3kb=round(float(phi[tss3 == 0].mean()), 5),
                         frac_tss3kb=round(float((tss3 == 1).mean()), 4), frac_t1na=round(float((t1 == 1).mean()), 4), phi_sd=round(float(phi.std()), 5), phi_cv=round(float(phi.std() / phi.mean()), 4),
                         phi_p99_p01=round(float(np.percentile(phi, 99) / np.percentile(phi, 1)), 2), phi_top1pct_mass=round(float(np.sort(phi)[-max(1, len(phi) // 100):].sum() / phi.sum()), 4)))
        print(rid, rows[-1]['rho_phi_dist'], rows[-1]['rho_phi_logdist'], rows[-1]['frac_dist_na'], flush=True)
import csv
with open(f'{OUT}/phi/phi_diag.tsv', 'w', newline='') as o:
    w = csv.DictWriter(o, fieldnames=list(rows[0]), delimiter='\t'); w.writeheader(); w.writerows(rows)
r_ = [r['rho_phi_dist'] for r in rows if not math.isnan(r['rho_phi_dist'])]
print('dist_tss in selected 18 cols:', 'dist_tss' in model_cols, '| regions', len(rows), 'rho median', round(float(np.median(r_)), 3), 'min', round(min(r_), 3), 'max', round(max(r_), 3), '| frac_dist_na median', round(float(np.median([r['frac_dist_na'] for r in rows])), 3))
