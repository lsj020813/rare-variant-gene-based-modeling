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
G=_config_path('${PROJECT_ROOT}/work/ref/groupfiles_chunks')
THR=2.5e-6
print("chrom  genes   TCHL_sig  TCHL_min      HTN_sig  any_p<1e-4")
rows=[]
for ch in range(1,23):
    exp=len(glob.glob(f'{G}/chr{ch}.part*.txt'))
    per={}
    ok=True
    for t in ('htn','dm','lip','tchl'):
        d=len(glob.glob(f'{R}/{t}.chr{ch}.part*.done'))
        if d!=exp or exp==0: ok=False
        ps=[]
        for f in glob.glob(f'{R}/{t}.chr{ch}.part[0-9][0-9][0-9]'):
            for line in open(f):
                v=line.rstrip("\n").split("\t")
                if len(v)<4 or v[0]=="Region": continue
                try: p=float(v[3])
                except ValueError: continue
                if 0<p<=1: ps.append(p)
        per[t]=ps
    if not per['tchl']: continue
    ng=len(per['tchl'])
    ts=sum(1 for p in per['tchl'] if p<THR); tm=min(per['tchl'])
    hs=sum(1 for p in per['htn'] if p<THR) if per['htn'] else 0
    a4=sum(sum(1 for p in v if p<1e-4) for v in per.values())
    rows.append((ch,ng,ts,tm,hs,a4,ok))
    print(f"chr{ch:<4}{ng:>6}{ts:>10}{tm:>12.1e}{hs:>13}{a4:>11}{'' if ok else '  (partial)'}")
print()
top=sorted(rows,key=lambda r:(-r[2],r[3]))[:3]
print("top by TCHL sig:", [(f"chr{r[0]}",r[2]) for r in top])
