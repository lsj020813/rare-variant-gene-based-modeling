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


import csv, json, sys, math
from pathlib import Path
sys.path.insert(0, _config_path('${PROJECT_ROOT}/work/run_l3'))
import l3_common as C

def load(path):
    d = {}
    with open(path, newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            d[r['Region']] = r
    return d

def main():
    arm, truth, out = sys.argv[1], sys.argv[2], Path(sys.argv[3])
    A, T = load(arm), load(truth)
    common = sorted(set(A) & set(T))
    cols = ['Pvalue', 'Pvalue_Burden', 'Pvalue_SKAT', 'BETA_Burden', 'SE_Burden', 'MAC',
            'Number_rare', 'Number_ultra_rare']
    exact = {c: 0 for c in cols}
    close = {c: 0 for c in cols}
    worst = {c: 0.0 for c in cols}
    for g in common:
        for c in cols:
            a, t = A[g].get(c, ''), T[g].get(c, '')
            if a == t:
                exact[c] += 1; close[c] += 1; continue
            try:
                fa, ft = float(a), float(t)
            except (TypeError, ValueError):
                continue
            rel = abs(fa - ft) / max(abs(ft), 1e-300)
            worst[c] = max(worst[c], rel)
            if rel <= 1e-9:
                close[c] += 1
    res = dict(arm_genes=len(A), truth_genes=len(T), common=len(common),
               only_in_arm=sorted(set(A) - set(T))[:20], only_in_truth=sorted(set(T) - set(A))[:20],
               identical_string=exact, identical_or_rel_1e9=close, worst_rel_diff=worst,
               gate_pass=bool(len(A) == len(T) == len(common) and close['Pvalue'] == len(common)))
    C.atomic_json(out / 'gate_none.json', res)
    print(json.dumps({k: res[k] for k in ('arm_genes', 'truth_genes', 'common', 'gate_pass')}))
    print('Pvalue identical:', exact['Pvalue'], '/ within 1e-9:', close['Pvalue'],
          '/ worst rel:', worst['Pvalue'])
    if not res['gate_pass']:
        print('GATE FAIL: G_none does not reproduce the truth file')
        raise SystemExit(1)
    print('GATE_NONE_PASS')

if __name__ == '__main__':
    main()
