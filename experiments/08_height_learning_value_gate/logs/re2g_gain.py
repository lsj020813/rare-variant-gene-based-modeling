
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import gzip,bisect,collections,glob
import pyarrow.parquet as pq, pyarrow.compute as pc
CH=("2","12")
G={c:[] for c in CH}; pc_ids=set()
for l in gzip.open(_config_path("${PROJECT_ROOT}/work/ref/deductive/gencode.nochr.gtf.gz"),"rt"):
    if l[0]=="#": continue
    a=l.split("\t")
    if a[2]!="gene" or 'gene_type "protein_coding"' not in a[8]: continue
    ch=a[0].replace("chr",""); gid=a[8].split('gene_id "')[1].split('"')[0].split(".")[0]; pc_ids.add(gid)
    if ch in CH: G[ch].append((int(a[3])-3000,int(a[4])+3000))
def merge(iv):
    iv=sorted(iv); m=[]
    for s,e in iv:
        if m and s<=m[-1][1]: m[-1][1]=max(m[-1][1],e)
        else: m.append([s,e])
    return m
B={c:merge(G[c]) for c in CH}
E={c:[] for c in CH}
for f in sorted(glob.glob(_config_path("${PROJECT_ROOT}/work/ref/re2g_all/ot_e2g/*.parquet"))):
    t=pq.read_table(f,columns=["chromosome","start","end","score","geneId"])
    t=t.filter(pc.and_(pc.is_in(pc.cast(t["chromosome"],"string"),value_set=__import__("pyarrow").array(["2","12","chr2","chr12"])),pc.greater_equal(t["score"],0.6)))
    for ch,s,e,g in zip(t["chromosome"].to_pylist(),t["start"].to_pylist(),t["end"].to_pylist(),t["geneId"].to_pylist()):
        if str(g).split(".")[0] in pc_ids: E[str(ch).replace("chr","")].append((int(s)+1,int(e)))
def inside(m,p):
    S=[x[0] for x in m]; i=bisect.bisect_right(S,p)-1
    return i>=0 and m[i][0]<=p<=m[i][1]
def bp(m): return sum(e-s+1 for s,e in m)
out={}
for c in CH:
    Em=merge(E[c]); U=merge(B[c]+Em)
    out[c]=dict(body_bp=bp(B[c]),union_bp=bp(U),added_bp=bp(U)-bp(B[c]),n_elem=len(E[c]))
    SB=[x[0] for x in B[c]]; SE=[x[0] for x in Em]
    o=0; rin=0
    for l in open(_config_path("${PROJECT_ROOT}/work/prs/out/hei/lead_snps.tsv")):
        a=l.rstrip().split("\t")
        if a[0]!=c: continue
        p=int(a[4])
        if not inside(B[c],p):
            o+=1; rin+=inside(Em,p)
    out[c].update(leads_outside_body=o, of_those_in_re2g=rin)
print(out)
