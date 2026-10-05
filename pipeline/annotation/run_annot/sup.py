
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import glob, math, json, os, re
D=_config_path("${PROJECT_ROOT}/work/ref/saige_step2_bwg"); OUT=_config_path("${PROJECT_ROOT}/work/ref/annot/supervision")
os.makedirs(OUT, exist_ok=True)
phs=("tchl","htn","dm","lip"); T={}; have={ph:set() for ph in phs}
for ph in phs:
    files=[f for f in glob.glob(f"{D}/{ph}.chr*") if f.split("/")[-1].count(".")==1]
    files+=[f for f in glob.glob(f"{D}_chr1/{ph}.chr1.part*") if re.search(r"part[0-9]{3}$", f)]
    for f in files:
        ch=f.split(".chr")[1].split(".")[0]; have[ph].add(ch)
        with open(f) as fh:
            hdr=next(fh).rstrip("\n").split("\t"); ix={h:i for i,h in enumerate(hdr)}
            for line in fh:
                a=line.rstrip("\n").split("\t")
                if len(a)!=len(hdr): continue
                try: p=float(a[ix["Pvalue"]]); pb=float(a[ix["Pvalue_Burden"]]); ps=float(a[ix["Pvalue_SKAT"]]); mac=float(a[ix["MAC"]])
                except: continue
                T.setdefault(a[0],{"chr":ch})[ph]=(p,pb,ps,int(mac))
missing={ph:sorted(set(map(str,range(1,23)))-have[ph],key=int) for ph in phs}
print("missing chroms:",missing)
genes=sorted(T); full=[g for g in genes if all(ph in T[g] for ph in phs)]
print("genes any:",len(genes),"| all 4:",len(full))
FOLD={"7":0,"13":0,"16":0,"19":0,"21":0,"2":1,"8":1,"9":1,"10":1,"15":1,"1":2,"6":2,"14":2,"3":3,"12":3,"17":3,"22":3,"4":4,"5":4,"11":4,"18":4,"20":4}
with open(f"{OUT}/gene_trait_skato.tsv","w") as o:
    o.write("gene\tgene_base\tchr\tfold\t"+"\t".join(f"{ph}_{k}" for ph in phs for k in ("p","pB","pS","mac"))+"\n")
    for g in full:
        o.write(f"{g}\t{g.split('.')[0]}\t{T[g]['chr']}\t{FOLD[T[g]['chr']]}\t"+"\t".join(f"{x:.6g}" for ph in phs for x in T[g][ph])+"\n")
def nl(p): return min(10.0,-math.log10(max(p,1e-300)))
X={ph:[nl(T[g][ph][0]) for g in full] for ph in phs}
def rank(a):
    idx=sorted(range(len(a)),key=lambda i:a[i]); r=[0.0]*len(a); i=0
    while i<len(a):
        j=i
        while j+1<len(a) and a[idx[j+1]]==a[idx[i]]: j+=1
        for k in range(i,j+1): r[idx[k]]=(i+j)/2
        i=j+1
    return r
def spear(a,b):
    RA,RB=rank(a),rank(b); ma=sum(RA)/len(a); mb=sum(RB)/len(b)
    num=sum((RA[i]-ma)*(RB[i]-mb) for i in range(len(a)))
    den=math.sqrt(sum((x-ma)**2 for x in RA)*sum((x-mb)**2 for x in RB)); return num/den
res={f"{a}-{b}":round(spear(X[a],X[b]),4) for i,a in enumerate(phs) for b in phs[i+1:]}
print("spearman(-log10 p, clip10):",json.dumps(res))
thr=0.05/18608
print("sig p<2.69e-6:",{ph:sum(1 for g in full if T[g][ph][0]<thr) for ph in phs})
print("frac -log10p>2:",{ph:round(sum(1 for x in X[ph] if x>2)/len(full),4) for ph in phs})
print("fold gene counts:",{f:sum(1 for g in full if FOLD[T[g]['chr']]==f) for f in range(5)})
json.dump({"n":len(full),"missing":missing,"spearman":res},open(f"{OUT}/cross_trait_corr.json","w"))
