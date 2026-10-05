import os
import json

def required(name):
    value = os.environ.get(name, "")
    if not value.strip():
        raise RuntimeError("Required environment variable is missing or empty: " + name)
    return value

PROJECT_ROOT = required("PROJECT_ROOT")
import glob, gzip, subprocess
G=f'{PROJECT_ROOT}/work/ref/groupfiles_chunks'
P=f'{PROJECT_ROOT}/work/ref/groupfiles_pilot'
BCF = os.environ.get("BCFTOOLS", "bcftools")
A=set()
for fp in glob.glob(f'{G}/chr19.part*.txt'):
    for line in open(fp):
        f=line.rstrip("\n").split()
        if len(f)>=3 and f[1]=="var": A.update(f[2:])
B=set()
for line in open(f'{P}/chr19.B_3kb_re2g.txt'):
    f=line.rstrip("\n").split()
    if len(f)>=3 and f[1]=="var": B.update(f[2:])
print("A distinct:", f"{len(A):,}", "| B distinct:", f"{len(B):,}")
print("A sample keys:", sorted(A)[:2])
print("B sample keys:", sorted(B)[:2])
print("A subset of B?", A <= B, "| B subset of A?", B <= A, "| intersection:", f"{len(A&B):,}")
txt = subprocess.run([BCF, "query", "-f", "%ID\n", f"{PROJECT_ROOT}/work/ref/lift38_keyed/chr19.keyed38.vcf.gz"],
                     check=True, capture_output=True, text=True).stdout
allk = {l[3:] if l.startswith("chr") else l for l in txt.split()}
print("keyed38 IDs total:", f"{len(allk):,}", "| A in that set:", f"{len(A&allk):,}", "| B in:", f"{len(B&allk):,}")
