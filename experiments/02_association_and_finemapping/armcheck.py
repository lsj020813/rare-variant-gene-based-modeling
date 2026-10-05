import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import gzip, glob, collections
GTF=_config_path('${PROJECT_ROOT}/work/ref/deductive/gencode.sorted.gtf.gz')
tss={}
with gzip.open(GTF,'rt') as fh:
    for line in fh:
        if line[0]=='#': continue
        f=line.split('\t')
        if f[2]!='gene' or f[0]!='chr19': continue
        if 'gene_type "protein_coding"' not in f[8]: continue
        i=f[8].find('gene_id "'); gid=f[8][i+9:f[8].find('"',i+9)].split('.')[0]
        tss[gid]= int(f[3]) if f[6]=='+' else int(f[4])
print("chr19 protein-coding genes:", len(tss))

dist=[]
for fp in glob.glob(_config_path('${PROJECT_ROOT}/work/ref/re2g/*.bed.gz')):
    with gzip.open(fp,'rt') as fh:
        for line in fh:
            f=line.rstrip('\n').split('\t')
            if len(f)<14 or f[0]!='chr19': continue
            g=f[13].split('.')[0]
            if g not in tss: continue
            mid=(int(f[1])+int(f[2]))//2
            dist.append(abs(mid-tss[g]))
dist.sort()
n=len(dist)
if n:
    import bisect
    for cut in (3000, 10000, 50000, 100000, 250000):
        k=bisect.bisect_right(dist,cut)
        print(f"  links within +/-{cut//1000:>3}kb: {k:>8,} / {n:,} = {k/n*100:5.1f}%")
    print(f"  median link distance: {dist[n//2]:,} bp | p90 {dist[int(n*0.9)]:,} bp | max {dist[-1]:,}")
