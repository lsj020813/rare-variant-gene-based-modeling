
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import csv, glob, os, collections
os.chdir(_config_path('${PROJECT_ROOT}/work/ref/saige_step2_bwg'))
files=[f for f in glob.glob('*') if os.path.isfile(f) and not f.endswith(('.done','.index','.log','_temp')) and 'singleAssoc' not in f]
THR=2.5e-6
omni=collections.defaultdict(set); skat=collections.defaultdict(set); burd=collections.defaultdict(set); n=collections.Counter()
for f in files:
    tr=f.split('.')[0]
    try: fh=open(f, encoding='utf-8', errors='replace')
    except OSError: continue
    for r in csv.DictReader(fh, delimiter='\t'):
        if r.get('Group')!='all': continue
        try:
            o=float(r['Pvalue']); s=float(r['Pvalue_SKAT']); b=float(r['Pvalue_Burden'])
        except (ValueError,KeyError,TypeError): continue
        g=r['Region']; n[tr]+=1
        if o<THR: omni[tr].add(g)
        if s<THR: skat[tr].add(g)
        if b<THR: burd[tr].add(g)
print('files', len(files), 'THR', THR)
print(f"{'trait':6s} {'tested':>9s} {'omni':>6s} {'skat':>6s} {'burden':>7s} {'skat_only':>10s} {'omni_only':>10s}")
TO=TS=TB=0
for tr in sorted(n):
    o,s,b=omni[tr],skat[tr],burd[tr]; TO+=len(o); TS+=len(s); TB+=len(b)
    print(f"{tr:6s} {n[tr]:9,} {len(o):6d} {len(s):6d} {len(b):7d} {len(s-o):10d} {len(o-s):10d}")
print(f"{'TOTAL':6s} {sum(n.values()):9,} {TO:6d} {TS:6d} {TB:7d}")
gain=[(tr,g) for tr in skat for g in sorted(skat[tr]-omni[tr])]
lost=[(tr,g) for tr in omni for g in sorted(omni[tr]-skat[tr])]
print('SKAT_ONLY_PAIRS', len(gain)); print(gain[:30])
print('OMNI_ONLY_PAIRS', len(lost)); print(lost[:15])
print('SKATCOUNT_DONE')
