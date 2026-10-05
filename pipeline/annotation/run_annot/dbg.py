
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import glob
D=_config_path("${PROJECT_ROOT}/work/ref/saige_step2_bwg")
for ph in ("tchl","htn","dm","lip"):
    fs=[f for f in glob.glob(f"{D}/{ph}.chr*") if f.split("/")[-1].count(".")==1]
    n=0; bad=0
    for f in fs:
        for line in open(f):
            a=line.rstrip("\n").split("\t")
            if a[0]=="Region" or len(a)<13: bad+=1; continue
            n+=1
    print(ph, "files", len(fs), "rows", n, "skipped", bad)
f=f"{D}/tchl.chr22"; a=open(f).readlines()[1].rstrip("\n").split("\t"); print("ncols tchl.chr22:", len(a), a[:3])
f=f"{D}/dm.chr22"; a=open(f).readlines()[1].rstrip("\n").split("\t"); print("ncols dm.chr22:", len(a), a[:3])
