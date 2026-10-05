import os
import json

def required(name):
    value = os.environ.get(name, "")
    if not value.strip():
        raise RuntimeError("Required environment variable is missing or empty: " + name)
    return value

PROJECT_ROOT = required("PROJECT_ROOT")
import glob
G=f'{PROJECT_ROOT}/work/ref/groupfiles_chunks'
pairs=0; uniq=set(); genes=set()
for fp in glob.glob(f'{G}/chr19.part*.txt'):
    for line in open(fp):
        f=line.rstrip("\n").split()
        if len(f)<3 or f[1]!="var": continue
        genes.add(f[0]); pairs+=len(f)-2; uniq.update(f[2:])
print(f"[A_100kb_re2g] genes {len(genes):,} | pairs {pairs:,} | distinct {len(uniq):,} | "
      f"mean m {pairs/len(genes):.1f} | redundancy {pairs/len(uniq):.2f}x")
