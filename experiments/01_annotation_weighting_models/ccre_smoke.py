
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import subprocess, json, collections
R=_config_path("${PROJECT_ROOT}/work/ref")
BT="bedtools"
CH="21"
out = subprocess.run(f"{BT} intersect -a {R}/b6_cards/chr{CH}.var.s.bed -b {R}/b6_cards/ccre.s.bed -wa -wb", 
                     shell=True, capture_output=True, text=True)
cls = collections.Counter()
seen = set()
for line in out.stdout.splitlines():
    f = line.split("\t")
    key = (f[0], f[1])
    cls[f[-1]] += 1
    seen.add(key)
nvar = int(subprocess.run(f"wc -l < {R}/b6_cards/chr{CH}.var.s.bed", shell=True, capture_output=True, text=True).stdout.strip())
lf = len(seen)/nvar
print(f"chr{CH}: variants {nvar:,}  labeled {len(seen):,} ({lf:.4f})")
print("classes:", dict(cls))
ref = json.load(open(f"{R}/b6_cards/chr{CH}.b6.json"))
print("b6 json labeled_frac:", ref["labeled_frac"], "| ccre_counts:", ref.get("ccre_counts"))
assert abs(lf - ref["labeled_frac"]) < 0.005, f"GATE FAIL: labeled frac {lf} vs {ref['labeled_frac']}"
print("SMOKE_CCRE_OK")
