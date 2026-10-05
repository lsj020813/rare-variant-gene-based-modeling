import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import glob, math
from statistics import median
R=_config_path('${PROJECT_ROOT}/work/ref/saige_step2_v2')
THR=2.5e-6
for t in ('htn','dm','lip','tchl'):
    ps=[]
    for f in glob.glob(f'{R}/{t}.chr*.part[0-9][0-9][0-9]'):
        for line in open(f):
            v=line.rstrip('\n').split('\t')
            if len(v)<4 or v[0]=='Region': continue
            try: p=float(v[3])
            except ValueError: continue
            if 0<p<=1: ps.append(p)
    if not ps:
        print(t, 'no rows'); continue
    lam=median([abs(math.log(p))*2 for p in ps])/1.3863
    print(f'{t}\tgenes={len(ps)}\tlambda={lam:.3f}\tsig={sum(1 for p in ps if p<THR)}\tsig1e4={sum(1 for p in ps if p<1e-4)}\tmin={min(ps):.2e}')
