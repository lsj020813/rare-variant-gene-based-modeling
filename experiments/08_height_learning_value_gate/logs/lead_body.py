
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import gzip,bisect,collections
G=collections.defaultdict(list)
for l in gzip.open(_config_path("${PROJECT_ROOT}/work/ref/deductive/gencode.nochr.gtf.gz"),"rt"):
    if l[0]=="#": continue
    a=l.split("\t")
    if a[2]!="gene" or 'gene_type "protein_coding"' not in a[8]: continue
    G[a[0].replace("chr","")].append((int(a[3])-3000,int(a[4])+3000))
M={}
for ch,iv in G.items():
    iv.sort(); m=[]
    for s,e in iv:
        if m and s<=m[-1][1]: m[-1][1]=max(m[-1][1],e)
        else: m.append([s,e])
    M[ch]=([x[0] for x in m],m)
def inside(ch,p):
    S,m=M[ch]; i=bisect.bisect_right(S,p)-1
    return i>=0 and m[i][0]<=p<=m[i][1]
c=collections.Counter()
for l in open(_config_path("${PROJECT_ROOT}/work/prs/out/hei/lead_snps.tsv")):
    a=l.rstrip().split("\t"); ch=a[0]; p=int(a[4])
    k="all"; c[(k,inside(ch,p))]+=1
    if ch in ("2","12"): c[("chr2_12",inside(ch,p))]+=1
print({k:v for k,v in sorted(c.items())})
for ch in ("2","12"):
    n=t=0
    for i,l in enumerate(gzip.open(_config_path(f'${{PROJECT_ROOT}}/work/fset/out/uni_chr{ch}.tsv.gz'),"rt")):
        if i==0: continue
        a=l.split("\t")
        if a[5]!="0" or float(a[3])<0.3: continue
        n+=1; t+=inside(ch,int(a[1]))
    print("chr",ch,"uni_noncoding",n,"in_body3kb",t,round(t/n,3))
