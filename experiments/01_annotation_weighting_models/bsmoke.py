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

import subprocess, collections, os, json
R=_config_path("${PROJECT_ROOT}/work/ref")
BCF="bcftools"
card={}
with open(f"{R}/l1/chr22.ccre_card.tsv") as fh:
    fh.readline()
    for line in fh:
        f=line.rstrip("\n").split("\t")
        card[f[1]]=f[2]
out=subprocess.run(f"{BCF} query -f '%POS\t%ID\n' {R}/lift38_keyed/chr22.keyed38.vcf.gz",
                   shell=True, capture_output=True, text=True)
k2p={}
for l in out.stdout.splitlines():
    p,kid=l.split("\t")
    k2p[kid[3:] if kid.startswith("chr") else kid]=str(int(p)-1)
genes=[]
for line in open(f"{R}/groupfiles_bwg/chr22.B_3kb_re2g.txt"):
    f=line.split()
    if f[1]=="var":
        genes.append((f[0], [k[3:] if k.startswith("chr") else k for k in f[2:]]))
    if len(genes)>=40: break
import numpy as np
samples=subprocess.run(f"{BCF} query -l {R}/band_vcf/chr22.band.vcf.gz", shell=True,
                       capture_output=True, text=True).stdout.split()
NS=len(samples)
assert NS==N_SAMPLES, f"GATE FAIL samples {NS}"
ok=0; part_fail=0
for g, keys in genes:
    pat=collections.defaultdict(list)
    for k in keys:
        p0=k2p.get(k)
        pat[card.get(p0,"none") if p0 else "none"].append(k)
    assert sum(len(v) for v in pat.values())==len(keys), f"GATE FAIL slots {g}"
    lo=min(int(k.split(":")[1]) for k in keys); hi=max(int(k.split(":")[1]) for k in keys)
    q=subprocess.run(f"{BCF} query -r 22:{lo}-{hi} -f '%CHROM:%POS:%REF:%ALT[\t%DS]\n' {R}/band_vcf/chr22.band.vcf.gz",
                     shell=True, capture_output=True, text=True)
    ds={}
    want=set(keys)
    for l in q.stdout.splitlines():
        f=l.split("\t")
        if f[0] in want:
            ds[f[0]]=np.array(f[1:], dtype=np.float32)
    assert len(ds)==len(keys), f"GATE FAIL DS fetch {g}: {len(ds)}/{len(keys)}"
    whole=np.zeros(NS, dtype=np.float32)
    psum=np.zeros(NS, dtype=np.float32)
    for k in keys: whole += ds[k]
    for kk, kl in pat.items():
        b=np.zeros(NS, dtype=np.float32)
        for k in kl: b += ds[k]
        psum += b
    if not np.allclose(whole, psum, atol=1e-3): part_fail += 1
    ok+=1
print(f"genes tested {ok}  partition-identity failures {part_fail}")
assert part_fail==0, "GATE FAIL: pattern partition identity"
print("B_SMOKE_OK")
