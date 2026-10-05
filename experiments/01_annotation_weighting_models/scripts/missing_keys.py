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


import glob, os, json
W=_config_path('${PROJECT_ROOT}/work/run_ourfm'); OUT=_config_path('${PROJECT_ROOT}/work/run_trackA')
regions={}
for t in ['dm','htn','lip','tchl']:
    for f in sorted(glob.glob(f'{W}/fm/{t}/*.ld.vars')):
        rid=os.path.basename(f)[:-len('.ld.vars')]
        if '.mi1000' in rid: continue
        regions[(t,rid)]=[l.strip() for l in open(f) if l.strip()]
assert len(regions)==81, len(regions)
def swap(k):
    a=k.split(':'); return f'{a[0]}:{a[1]}:{a[3]}:{a[2]}'
have=set()
for f in glob.glob(f'{OUT}/annot/*.bbj.tsv'):
    with open(f) as fh:
        hdr=fh.readline().rstrip('\n').split('\t'); assert hdr[0]=='orient' and hdr[1]=='variant_hg19'
        for line in fh:
            a=line.split('\t',2); k=a[1]
            have.add(k if a[0]=='0' else swap(k))
missing=set(); per={}
for (t,rid),v in regions.items():
    m=[k for k in v if k not in have]; missing.update(m); per[rid]=dict(n_fm=len(v),n_missing=len(m))
os.makedirs(f'{OUT}/annot_trA',exist_ok=True)
with open(f'{OUT}/annot_trA/missing_keys37.txt','w') as o:
    for k in sorted(missing): o.write(k+'\n')
tot=sum(p['n_fm'] for p in per.values()); mis=sum(p['n_missing'] for p in per.values())
json.dump(dict(n_regions=len(regions),n_var_total=tot,n_missing_total=mis,n_missing_unique=len(missing),have_keys=len(have),per_region=per),open(f'{OUT}/annot_trA/missing_keys.json','w'),indent=1)
print('n_regions',len(regions),'n_var',tot,'have_keys',len(have),'missing_total',mis,'missing_unique',len(missing),'frac',round(mis/tot,4))
open(f'{OUT}/annot_trA/missing_keys37.txt.done','w').write(f'ok {len(missing)}\n')
