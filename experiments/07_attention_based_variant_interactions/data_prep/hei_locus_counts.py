import os as _os
if not _os.environ.get("PROJECT_ROOT", "").strip():
    raise ValueError("PROJECT_ROOT must be set and nonblank")

import gzip, zipfile, io, json
Z=(_os.environ["PROJECT_ROOT"] + '/work/ref/bbj/hum0197.v3.BBJ.Hei.v1.zip')
f=io.TextIOWrapper(gzip.GzipFile(fileobj=zipfile.ZipFile(Z).open("hum0197.v3.BBJ.Hei.v1/GWASsummary_Height_Japanese_SakaueKanai2020.auto.txt.gz")))
h=f.readline().rstrip("\n").split("\t"); ix={k:i for i,k in enumerate(h)}
sig=[]; n=0; chi=[]
for line in f:
    p=line.split("\t"); n+=1; pv=float(p[ix["P_BOLT_LMM_INF"]]); ch=p[ix["CHR"]]; bp=int(p[ix["BP"]])
    if n%200==0: chi.append(pv)
    if pv<5e-8 and not (ch=="6" and 25_000_000<=bp<=34_000_000): sig.append((pv,ch,bp))
sig.sort(); res={"n_records":n,"n_sig":len(sig)}
for sp in (250000,500000):
    ch_=[]
    for pv,ch,bp in sig:
        if all(not (c==ch and abs(b-bp)<=sp) for _,c,b in ch_): ch_.append((pv,ch,bp))
    wins=sorted(((ch,max(0,bp-250000),bp+250000) for _,ch,bp in ch_), key=lambda x:(int(x[0]),x[1])); m=[]
    for ch,s,e in wins:
        if m and m[-1][0]==ch and s<=m[-1][2]: m[-1][2]=max(m[-1][2],e)
        else: m.append([ch,s,e])
    res[f"leads_{sp//1000}kb"]=len(ch_); res[f"loci_merged_{sp//1000}kb"]=len(m)
chi.sort(); res["median_p_sample"]=chi[len(chi)//2]
json.dump(res,open((_os.environ["PROJECT_ROOT"] + '/work/prs/out/hei/locus_rule_counts.json'),"w")); print(res)
