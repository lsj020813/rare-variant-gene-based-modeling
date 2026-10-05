import os
import json

def required(name):
    value = os.environ.get(name, "")
    if not value.strip():
        raise RuntimeError("Required environment variable is missing or empty: " + name)
    return value

PROJECT_ROOT = required("PROJECT_ROOT")
import gzip, glob
from collections import defaultdict
elem = defaultdict(set)
for fp in sorted(glob.glob(f'{PROJECT_ROOT}/work/ref/re2g/*.bed.gz')):
    with gzip.open(fp,'rt') as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 14: continue
            elem[(f[0], f[1], f[2])].add(f[13].split('.')[0])
n = len(elem)
dist = defaultdict(int)
for k, gs in elem.items(): dist[len(gs)] += 1
print(f"distinct rE2G elements (10-tissue union): {n:,}")
multi = sum(v for k, v in dist.items() if k >= 2)
print(f"  targeting 1 gene : {dist[1]:,} ({dist[1]/n*100:.1f}%)")
print(f"  targeting >=2    : {multi:,} ({multi/n*100:.1f}%)")
for k in sorted(dist):
    if k <= 8 or k == max(dist): print(f"    {k:>2} genes: {dist[k]:,}")
print(f"  max genes for one element: {max(dist)}")
tot_links = sum(k*v for k,v in dist.items())
print(f"  total element-gene links: {tot_links:,} | mean genes/element {tot_links/n:.2f}")
