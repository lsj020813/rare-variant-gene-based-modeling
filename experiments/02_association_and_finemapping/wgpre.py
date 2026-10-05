
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import subprocess, glob, os, gzip
BCF="bcftools"
R=_config_path("${PROJECT_ROOT}/work/ref")
checks=[]
def ck(n,v,ok): checks.append((n,str(v)[:40],ok))

print("### WG B-arm inputs ###")
ck("band VCF chr1..22", sum(os.path.exists(f"{R}/band_vcf/chr{i}.band.vcf.gz") for i in range(1,23)), 
   all(os.path.exists(f"{R}/band_vcf/chr{i}.band.vcf.gz") for i in range(1,23)))
b_pilot = glob.glob(f"{R}/groupfiles_pilot/chr*.B_3kb_re2g.txt")
ck("B-arm group files present", f"{len(b_pilot)} files: {[os.path.basename(x).split('.')[0] for x in b_pilot]}", len(b_pilot)>0)
for t in ("htn","dm","lip","tchl"):
    ok=os.path.exists(f"{R}/saige_step1_v4/{t}_v4.rda") and os.path.exists(f"{R}/saige_step1_v4/{t}_v4.varianceRatio.txt")
    ck(f"{t} null+VR (v4)", "present" if ok else "MISSING", ok)
ck("pilot builder pg.py", os.path.exists(_config_path("${PROJECT_ROOT}/work/run_pilot/pg.py")), os.path.exists(_config_path("${PROJECT_ROOT}/work/run_pilot/pg.py")))
if os.path.exists(_config_path("${PROJECT_ROOT}/work/run_pilot/pg.py")):
    src=open(_config_path("${PROJECT_ROOT}/work/run_pilot/pg.py")).read()
    ck("builder CH hardcoded", "CH=19" if 'CH = "19"' in src else "param", 'CH = "19"' in src)
    ck("builder anno label", "all" if '["all"]' in src else "OTHER", '["all"]' in src)
    ck("builder assign rule B", "3kb+rE2G" if "window_assign(3_000, True)" in src else "?", "window_assign(3_000, True)" in src)

n_re2g = len(glob.glob(f"{R}/re2g/*.bed.gz"))
ck("rE2G tissue files", n_re2g, n_re2g>0)
n_keyed = len(glob.glob(f"{R}/lift38_keyed/chr*.keyed38.vcf.gz"))
ck("keyed38 maps chr1..22", n_keyed, n_keyed>=22)

print(f"{'check':<32}{'value':<42}ok")
for n,v,ok in checks: print(f"{n:<32}{v:<42}{'OK' if ok else '** FAIL **'}")
print("FAILS:", sum(1 for _,_,ok in checks if not ok))
print("WG_PREFLIGHT_DONE")
