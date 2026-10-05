import os as _cfg_os
import re as _cfg_re
import math as _cfg_math

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)

def _config_number(name, cast=int, positive=False):
    value = cast(_cfg_os.environ[name])
    if not _cfg_math.isfinite(value) or (positive and value <= 0):
        raise ValueError(name + " must be finite" + (" and positive" if positive else ""))
    return value

import os as _os
N_SAMPLES = _config_number("N_SAMPLES", int, True)

import subprocess, glob, os, collections
BCF="bcftools"
R=_config_path("${PROJECT_ROOT}/work/ref")
RAW=_config_path("${GENOTYPE_DIR}")

checks=[]
def ck(name, val, ok): checks.append((name, val, ok)); return ok

for ch in ("1","5","16","19"):
    f=f"{RAW}/chr{ch}.vcf.gz"
    ctg=subprocess.run(f"{BCF} query -f '%CHROM\n' -r {ch}:1-2000000 {f} 2>/dev/null | head -1",
                       shell=True,capture_output=True,text=True).stdout.strip()
    if not ctg:
        ctg=subprocess.run(f"zcat {f} | grep -v '^#' | head -1 | cut -f1", shell=True,capture_output=True,text=True).stdout.strip()
    ck(f"raw chr{ch} contig", ctg, ctg==ch)

f=f"{RAW}/chr19.vcf.gz"
hdr=subprocess.run(f"{BCF} view -h {f} | grep '##FORMAT'", shell=True,capture_output=True,text=True).stdout
ck("raw has FORMAT/DS", "ID=DS" in hdr, "ID=DS" in hdr)
ck("raw has FORMAT/GT", "ID=GT" in hdr, "ID=GT" in hdr)

ids=subprocess.run(f"{BCF} query -r 19:44900000-44910000 -f '%ID\n' {f}", shell=True,capture_output=True,text=True).stdout.split()
ck("raw ID field", (ids[:2] if ids else "none"), all(i=="." for i in ids[:20]) if ids else False)

bf=f"{R}/band_vcf/chr19.band.vcf.gz"
bctg=subprocess.run(f"zcat {bf} | grep -v '^#' | head -1 | cut -f1", shell=True,capture_output=True,text=True).stdout.strip()
ck("band chr19 contig", bctg, bctg=="19")
ns_raw=subprocess.run(f"{BCF} query -l {f} | wc -l", shell=True,capture_output=True,text=True).stdout.strip()
ns_band=subprocess.run(f"{BCF} query -l {bf} | wc -l", shell=True,capture_output=True,text=True).stdout.strip()
ck("raw samples", ns_raw, ns_raw==str(N_SAMPLES))
ck("band samples", ns_band, ns_band==str(N_SAMPLES))
s_raw=subprocess.run(f"{BCF} query -l {f}", shell=True,capture_output=True,text=True).stdout.split()
s_band=subprocess.run(f"{BCF} query -l {bf}", shell=True,capture_output=True,text=True).stdout.split()
ck("sample order identical", f"nonempty={bool(s_raw) and bool(s_band)}", bool(s_raw) and s_raw==s_band)

for t in ("htn","dm","lip","tchl"):
    ok = os.path.getsize(f"{R}/saige_step1_v4/{t}_v4.rda")>0 and os.path.getsize(f"{R}/saige_step1_v4/{t}_v4.varianceRatio.txt")>0
    ck(f"{t} null model + VR", "present", ok)

for ch in ("1","5","16","19"):
    n=len(glob.glob(f"{R}/groupfiles_chunks/chr{ch}.part*.txt"))
    ck(f"groupfiles chr{ch}", f"{n} chunks", n>0)

rec=subprocess.run(f"{BCF} query -r 19:45411941-45411941 -f '%CHROM:%POS_%REF/%ALT\n' {f}",
                   shell=True,capture_output=True,text=True).stdout.strip()
ck("example condition id (rs429358)", rec, bool(rec))

print(f"{'check':<34}{'value':<28}ok")
bad=0
for n,v,ok in checks:
    print(f"{n:<34}{str(v)[:27]:<28}{'OK' if ok else '** FAIL **'}")
    if not ok: bad+=1
print(f"\nFAILS: {bad}")
print("PREFLIGHT_DONE")
