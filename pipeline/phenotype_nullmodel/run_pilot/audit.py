import os
import json

def required(name):
    value = os.environ.get(name, "")
    if not value.strip():
        raise RuntimeError("Required environment variable is missing or empty: " + name)
    return value

PROJECT_ROOT = required("PROJECT_ROOT")
import glob, os, re
R=f'{PROJECT_ROOT}/work/ref/saige_step2_v2'
G=f'{PROJECT_ROOT}/work/ref/groupfiles_chunks'
traits=('htn','dm','lip','tchl')
missing=[]; empty=[]; nodone=[]; ok=0; genes_by_trait={t:0 for t in traits}
for ch in range(1,23):
    parts=sorted(os.path.basename(p) for p in glob.glob(f'{G}/chr{ch}.part*.txt'))
    for t in traits:
        for p in parts:
            tag=p.replace('.txt','').replace('chr','chr')
            f=f'{R}/{t}.{tag}'
            if not os.path.exists(f): missing.append(f'{t}.{tag}'); continue
            n=sum(1 for l in open(f) if not l.startswith('Region'))
            if n==0: empty.append(f'{t}.{tag}')
            elif not os.path.exists(f+'.done'): nodone.append(f'{t}.{tag}'); genes_by_trait[t]+=n
            else: ok+=1; genes_by_trait[t]+=n
exp_parts=len(glob.glob(f'{G}/chr*.part*.txt'))
print(f"expected chunk-runs: {exp_parts*4} ({exp_parts} chunks x 4 traits)")
print(f"  complete (result+done): {ok}")
print(f"  result present, no .done marker: {len(nodone)} {nodone[:6]}")
print(f"  EMPTY result: {len(empty)} {empty[:6]}")
print(f"  MISSING entirely: {len(missing)} {missing[:8]}")
print("\ntested genes by trait:", genes_by_trait)
print("total:", sum(genes_by_trait.values()))
