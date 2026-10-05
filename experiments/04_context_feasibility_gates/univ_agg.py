
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import glob, collections, json, csv
rows=[]
for f in glob.glob(_config_path("${PROJECT_ROOT}/work/gate1/out/universe_*_chr*.tsv")):
    for line in open(f):
        p=line.rstrip("\n").split("\t")
        if len(p)==7 and p[2]!="TOTAL": rows.append(p)
agg=collections.Counter(); tot=collections.Counter()
for tag,ch,b,rb,t,vt,al,*rest in [r+[None] for r in rows]:
    n=int(rows[0][-1]) if False else None
for p in rows:
    tag,ch,b,rb,t,vt,al = p[0],p[1],p[2],p[3],p[4],p[5],p[6]
for p in rows: pass
agg=collections.Counter(); 
for line_f in glob.glob(_config_path("${PROJECT_ROOT}/work/gate1/out/universe_*_chr*.tsv")):
    for line in open(line_f):
        p=line.rstrip("\n").split("\t")
        if len(p)!=8: continue
        tag,ch,b,rb,t,vt,al,n = p
        if b=="TOTAL": continue
        agg[(tag,b,rb,t,vt,al)] += int(n)
with open(_config_path("${PROJECT_ROOT}/work/gate1/out/tables/universe_qc_census.csv"),"w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["source","maf_bin","r2_band","typed","var_type","allelic","n"])
    for k,v in sorted(agg.items()): w.writerow(list(k)+[v])
byb=collections.Counter(); 
for (tag,b,rb,t,vt,al),v in agg.items(): byb[(tag,b)]+=v
print(json.dumps({f"{k[0]}|{k[1]}":v for k,v in sorted(byb.items())}, indent=0))
