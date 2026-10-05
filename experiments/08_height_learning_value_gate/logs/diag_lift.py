
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import gzip
UNI=_config_path("${PROJECT_ROOT}/work/fset/out/uni_chr12.tsv.gz")
GTF=_config_path("${PROJECT_ROOT}/work/ref/deductive/gencode.nochr.gtf.gz")
def chrom(v):
    v=str(v).strip(); v=v[3:] if v.lower().startswith("chr") else v
    return str(int(v)) if v.isdigit() else v
genes=[]
for l in gzip.open(GTF,"rt"):
    if l[0]=="#": continue
    a=l.split("\t")
    if a[2]!="gene" or 'gene_type "protein_coding"' not in a[8]: continue
    if a[0].replace("chr","")!="12": continue
    genes.append((max(1,int(a[3])-3000),int(a[4])+3000))
genes.sort()
def inside(p):
    import bisect
    S=[g[0] for g in genes]; i=bisect.bisect_right(S,p)-1
    return i>=0 and genes[i][0]<=p<=genes[i][1]
keys={}
with gzip.open(UNI,"rt") as f:
    h=f.readline().rstrip("\n").split("\t"); ix={c:i for i,c in enumerate(h)}
    for l in f:
        a=l.rstrip("\n").split("\t")
        if a[ix["is_coding"]]!="0" or float(a[ix["r2"]])<0.3: continue
        k=a[ix["key"]]; c,p,r,al=k.split(":")
        p=int(p)
        if not inside(p): continue
        keys[(chrom(c),p,r.upper(),al.upper())]=1
print("assigned keys", len(keys))
import re
def rc(s): return s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]
def orient(sr,sa,tr,ta):
    for r,a,sign,lab in ((sr,sa,1,"direct"),(sa,sr,-1,"swap"),(rc(sr),rc(sa),1,"rc"),(rc(sa),rc(sr),-1,"swap_rc")):
        if (r,a)==(tr,ta): return sign,lab
    return None
bad=[]; nmulti=0
with gzip.open(_config_path("${PROJECT_ROOT}/work/ref/lift38_keyed/chr12.keyed38.vcf.gz"),"rt") as h:
    for line in h:
        if line[0]=="#": continue
        row=line.rstrip("\n").split("\t",8)
        fields=row[2].split(":")
        if len(fields)!=4 or not fields[1].isdigit(): continue
        k=(chrom(fields[0]),int(fields[1]),fields[2].upper(),fields[3].upper())
        if k not in keys: continue
        alts=row[4].upper().split(",")
        if len(alts)>1: nmulti+=1
        matches=[m for m in (orient(k[2],k[3],row[3].upper(),a) for a in alts) if m is not None]
        if len(matches)!=1: bad.append((k,row[3],alts,matches))
print("multi_alt_rows_among_wanted", nmulti, "ambiguous_or_unmatched", len(bad))
for x in bad[:5]: print(x)
