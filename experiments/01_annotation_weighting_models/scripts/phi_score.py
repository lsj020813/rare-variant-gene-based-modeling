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


import sys, os, glob, gzip, json, time, math
sys.path.insert(0, _config_path('${PROJECT_ROOT}/work/run_band15/model_v10_out'))
import l1_train_v10 as L
L.load_libraries(4)
import numpy as np
W = _config_path('${PROJECT_ROOT}/work/run_ourfm'); OUT = _config_path('${PROJECT_ROOT}/work/run_trackA')
MODEL = _config_path('${PROJECT_ROOT}/work/run_band15/model_v10_out/e6_runs/primary_f0/phi.json')
only = sys.argv[1].split(',') if len(sys.argv) > 1 and sys.argv[1] else None
allow_no_trA = '--allow-no-trA' in sys.argv
model = json.load(open(MODEL)); cols = model['design']['cols']
pred = L.PhiPredictor(model)
os.makedirs(f'{OUT}/phi', exist_ok=True)
t0 = time.time()
def num(v):
    return float('nan') if v in ('', 'NA', 'nan') else float(v)
trA = {}
fp = f'{OUT}/annot_trA/trA_annot_missing.tsv.gz'
if os.path.exists(fp + '.done'):
    with gzip.open(fp, 'rt') as fh:
        hdr = fh.readline().rstrip('\n').split('\t'); ci = {c: hdr.index(c) for c in cols}; ki = hdr.index('variant_hg19')
        for line in fh:
            a = line.rstrip('\n').split('\t'); trA[a[ki]] = [num(a[ci[c]]) for c in cols]
    print('trA_rows', len(trA))
else:
    assert allow_no_trA, 'GATE FAIL: trA annot not done (use --allow-no-trA for smoke)'
    print('trA_rows 0 (smoke: no trA)')
def swap(k):
    a = k.split(':'); return f'{a[0]}:{a[1]}:{a[3]}:{a[2]}'
rows = []
for t in ['dm', 'htn', 'lip', 'tchl']:
    for f in sorted(glob.glob(f'{W}/fm/{t}/*.ld.vars')):
        rid = os.path.basename(f)[:-len('.ld.vars')]
        if '.mi1000' in rid: continue
        if only and rid not in only: continue
        keys = [l.strip() for l in open(f) if l.strip()]
        bbj = {}
        bf = f'{OUT}/annot/{rid}.bbj.tsv'
        if os.path.exists(bf):
            with open(bf) as fh:
                hdr = fh.readline().rstrip('\n').split('\t'); ci = {c: hdr.index(c) for c in cols}; ki = hdr.index('variant_hg19')
                for line in fh:
                    a = line.rstrip('\n').split('\t'); k = a[ki] if a[0] == '0' else swap(a[ki])
                    bbj[k] = [num(a[ci[c]]) for c in cols]
        X = []; src = []
        for k in keys:
            if k in bbj: X.append(bbj[k]); src.append('bbj')
            elif k in trA: X.append(trA[k]); src.append('trA')
            else: X.append([float('nan')] * len(cols)); src.append('median')
        X = np.asarray(X, dtype=float); src = np.asarray(src)
        have = src != 'median'
        phi = np.full(len(keys), np.nan)
        if have.sum() == 0: raise RuntimeError(f'GATE FAIL: {rid} 주석 있는 변이 0')
        phi[have] = pred.predict(X[have], cols)
        med = float(np.median(phi[have])); phi[~have] = med
        assert np.isfinite(phi).all() and (phi > 0).all()
        t1na = X[have, cols.index('t1_na')] == 1
        with open(f'{OUT}/phi/{rid}.phi.tsv', 'w') as o:
            o.write('variant_hg19\tphi\tsrc\n')
            for k, p, s_ in zip(keys, phi, src): o.write(f'{k}\t{p:.6g}\t{s_}\n')
        r = dict(trait=t, region=rid, n_fm=len(keys), n_bbj=int((src == 'bbj').sum()), n_trA=int((src == 'trA').sum()), n_median=int((src == 'median').sum()),
                 frac_median=round(float((src == 'median').mean()), 4), frac_t1na_of_annotated=round(float(t1na.mean()), 4), phi_median=round(med, 5),
                 phi_min=round(float(phi.min()), 5), phi_max=round(float(phi.max()), 5), phi_mean=round(float(phi.mean()), 5),
                 phi_ratio_max_min=round(float(phi.max() / phi.min()), 2), phi_cv=round(float(phi.std() / phi.mean()), 4))
        rows.append(r)
        print(rid, r['n_fm'], r['n_bbj'], r['n_trA'], r['n_median'], r['frac_median'], r['phi_median'], r['phi_ratio_max_min'], flush=True)
import csv
tag = 'smoke' if only else 'all'
with open(f'{OUT}/phi/phi_coverage_{tag}.tsv', 'w', newline='') as o:
    w = csv.DictWriter(o, fieldnames=list(rows[0]), delimiter='\t'); w.writeheader(); w.writerows(rows)
tot = sum(r['n_fm'] for r in rows); nm = sum(r['n_median'] for r in rows); nt = sum(r['n_trA'] for r in rows); nb = sum(r['n_bbj'] for r in rows)
print('TOTAL regions', len(rows), 'n_fm', tot, 'bbj', nb, 'trA', nt, 'median', nm, 'frac_median', round(nm / tot, 4), 'elapsed', round(time.time() - t0))
if not only: open(f'{OUT}/phi/PHI_DONE', 'w').write(f'ok {len(rows)} {tot} {nm}\n')
