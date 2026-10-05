#!/usr/bin/env python
import os as _os
if not _os.environ.get("PROJECT_ROOT", "").strip():
    raise ValueError("PROJECT_ROOT must be set and nonblank")
import gzip, zipfile, io, json, os, numpy as np
os.makedirs((_os.environ["PROJECT_ROOT"] + '/work/prs/out/ht'), exist_ok=True)
Z=(_os.environ["PROJECT_ROOT"] + '/work/ref/bbj/hum0197.v3.BBJ.Ht.v1.zip')
W=(_os.environ["PROJECT_ROOT"] + '/work/prs')
zf=zipfile.ZipFile(Z); f=io.TextIOWrapper(gzip.GzipFile(fileobj=zf.open("hum0197.v3.BBJ.Ht.v1/GWASsummary_Ht_Japanese_SakaueKanai2020.auto.txt.gz")))
h=f.readline().rstrip("\n").split("\t"); ix={k:i for i,k in enumerate(h)}
out=open(W+"/out/bbj_ht_prscs.txt.tmp","w"); out.write("SNP\tA1\tA2\tBETA\tP\n")
sig=[]; nest=[]; n=0
for line in f:
    p=line.rstrip("\n").split("\t")
    rs=p[ix["SNP"]]; ch=p[ix["CHR"]]; bp=int(p[ix["BP"]])
    a1=p[ix["ALLELE1"]]; a0=p[ix["ALLELE0"]]; fr=float(p[ix["A1FREQ"]])
    b=float(p[ix["BETA"]]); se=float(p[ix["SE"]]); pv=float(p[ix["P_BOLT_LMM_INF"]])
    n+=1
    if rs.startswith("rs"): out.write(f"{rs}\t{a1}\t{a0}\t{b}\t{pv}\n")
    if 0.05<fr<0.95 and se>0 and n%50==0: nest.append(1/(2*fr*(1-fr)*se*se))
    if pv<5e-8 and not (ch=="6" and 25_000_000<=bp<=34_000_000): sig.append((pv,ch,bp))
out.close(); import os; os.replace(W+"/out/bbj_ht_prscs.txt.tmp", W+"/out/bbj_ht_prscs.txt")
sig.sort(); chosen=[]
for pv,ch,bp in sig:
    if all(not (c==ch and abs(b-bp)<=500000) for _,c,b in chosen): chosen.append((pv,ch,bp))
n_leads_all=len(chosen); chosen=chosen[:42]
wins=sorted((ch,max(0,bp-250000),bp+250000) for _,ch,bp in chosen)
merged=[]
for ch,s,e in sorted(wins, key=lambda x:(int(x[0]),x[1])):
    if merged and merged[-1][0]==ch and s<=merged[-1][2]: merged[-1][2]=max(merged[-1][2],e)
    else: merged.append([ch,s,e])
with open(W+"/out/ht/loci.bed","w") as g:
    for i,(ch,s,e) in enumerate(merged): g.write(f"{ch}\t{s}\t{e}\tL{i:03d}\n")
lead=sorted((ch,bp) for _,ch,bp in chosen)
with open(W+"/out/ht/lead_snps.tsv","w") as g:
    for ch,bp in lead: g.write(f"{ch}\t{bp}\n")
res=dict(n_leads_all_500kb=n_leads_all, n_records=n, n_sig=len(sig), n_leads=len(chosen), n_loci=len(merged), N_est_median=float(np.median(nest)), N_est_iqr=[float(np.percentile(nest,25)),float(np.percentile(nest,75))],
         loci_kb_median=float(np.median([(e-s)/1e3 for _,s,e in merged])))
json.dump(res, open(W+"/out/ht/bbj_prep_summary.json","w"), indent=1); print(json.dumps(res))
