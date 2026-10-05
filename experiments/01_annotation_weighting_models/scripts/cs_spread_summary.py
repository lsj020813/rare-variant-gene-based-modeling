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


import csv, json, math, statistics as stt, os
import numpy as np
W = _config_path('${PROJECT_ROOT}/work/run_trackA')
rows = list(csv.DictReader(open(f'{W}/trackA_cs_phi_spread.tsv'), delimiter='\t'))
met = {r['region']: r for r in csv.DictReader(open(f'{W}/trackA_region_metrics_v2.tsv'), delimiter='\t')}
def f(x): return float('nan') if x in ('NA', '') else float(x)
def q(v): v = [x for x in v if not math.isnan(x)]; return dict(n=len(v), median=float(np.median(v)), q1=float(np.percentile(v, 25)), q3=float(np.percentile(v, 75)), min=min(v), max=max(v)) if v else dict(n=0)
multi = [r for r in rows if int(r['cs_size']) >= 2]
def block(rs):
    return dict(n_cs=len(rs), ratio_max_min=q([f(r['cs_phi_ratio_max_min']) for r in rs]), p99_p01=q([f(r['cs_phi_p99_p01']) for r in rs]), cv=q([f(r['cs_phi_cv']) for r in rs]),
                ratio_within_over_region=q([f(r['cs_phi_ratio_max_min']) / f(r['region_phi_ratio_max_min']) for r in rs]), cs_size=q([float(r['cs_size']) for r in rs]), min_abs_corr=q([f(r['min_abs_corr']) for r in rs]))
S = dict(n_rows=len(rows), n_regions=len({r['region']: 1 for r in rows}), n_regions_alpha=len({r['region'] for r in rows if r['src'] == 'alpha'}), n_cs_size1=sum(1 for r in rows if r['cs_size'] == '1'),
         all_multi=block(multi), pure_multi=block([r for r in multi if r['pure'] == 'TRUE']), impure_multi=block([r for r in multi if r['pure'] == 'FALSE']),
         by_size_bin={b: block([r for r in multi if r['size_bin'] == b]) for b in ['2-4', '5-9', '10-29', '30+']},
         by_size_bin_pure={b: block([r for r in multi if r['size_bin'] == b and r['pure'] == 'TRUE']) for b in ['2-4', '5-9', '10-29', '30+']})
shr = [r for r in multi if r['pure'] == 'TRUE' and met.get(r['region'], {}).get('A_pct_cs_shrink', 'NA') != 'NA' and f(met[r['region']]['A_pct_cs_shrink']) >= .95]
nsh = [r for r in multi if r['pure'] == 'TRUE' and met.get(r['region'], {}).get('A_pct_cs_shrink', 'NA') != 'NA' and f(met[r['region']]['A_pct_cs_shrink']) < .95]
S['cs_shrink_regions_pure_cs'] = dict(regions=sorted({r['region'] for r in shr}), **block(shr)) if shr else dict(n_cs=0)
S['no_shrink_regions_pure_cs'] = block(nsh)
exc = [r for r in multi if r['pure'] == 'TRUE' and met.get(r['region'], {}).get('A_exceeds', 'NA') == '1']; nex = [r for r in multi if r['pure'] == 'TRUE' and met.get(r['region'], {}).get('A_exceeds', 'NA') == '0']
S['A_exceeds_regions_pure_cs'] = block(exc); S['A_not_exceeds_regions_pure_cs'] = block(nex)
from scipy.stats import spearmanr
x = [f(r['min_abs_corr']) for r in multi]; y = [f(r['cs_phi_ratio_max_min']) for r in multi]; z = [math.log(float(r['cs_size'])) for r in multi]
S['spearman_minabscorr_vs_ratio'] = float(spearmanr(x, y).correlation); S['spearman_logsize_vs_ratio'] = float(spearmanr(z, y).correlation)
med = S['pure_multi']['ratio_max_min'].get('median', float('nan')); med_all = S['all_multi']['ratio_max_min'].get('median', float('nan'))
rule = lambda m: '<2.0: 사전으로 CS 경계를 넘는 것이 산술적으로 거의 불가능(기계적 원인 확정)' if m < 2 else ('2.0~5.0: 가능하나 약함 — 원인은 φ 값의 정확성 쪽' if m <= 5 else '>5.0: 산술적 여지 충분 — φ 가 틀린 변이를 가리킴(더 강한 음성)')
S['verdict_rule'] = dict(pure_cs_median_ratio=med, pure=rule(med), all_cs_median_ratio=med_all, all=rule(med_all))
json.dump(S, open(f'{W}/trackA_cs_phi_spread.summary.json', 'w'), indent=1, ensure_ascii=False)
p = f'{W}/trackA_summary_v2.json'; s2 = json.load(open(p)); s2['cs_phi_spread'] = S; json.dump(s2, open(p, 'w'), indent=1, ensure_ascii=False)
print(json.dumps(dict(n_rows=S['n_rows'], n_regions=S['n_regions'], size1=S['n_cs_size1'], pure=S['pure_multi']['ratio_max_min'], impure=S['impure_multi']['ratio_max_min'], within_over_region=S['pure_multi']['ratio_within_over_region'], by_bin={b: S['by_size_bin'][b]['ratio_max_min'].get('median') for b in S['by_size_bin']}, rho_corr=S['spearman_minabscorr_vs_ratio'], rho_size=S['spearman_logsize_vs_ratio'], verdict=S['verdict_rule']), ensure_ascii=False))
