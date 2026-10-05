#!/usr/bin/env python
import os as _os
if not _os.environ.get("PROJECT_ROOT", "").strip():
    raise ValueError("PROJECT_ROOT must be set and nonblank")
import gzip, zipfile, io, json, numpy as np
Z=(_os.environ["PROJECT_ROOT"] + '/work/ref/bbj/hum0197.v3.BBJ.Hei.v1.zip')
W=(_os.environ["PROJECT_ROOT"] + '/work/prs')
zf=zipfile.ZipFile(Z); f=io.TextIOWrapper(gzip.GzipFile(fileobj=zf.open("hum0197.v3.BBJ.Hei.v1/GWASsummary_Height_Japanese_SakaueKanai2020.auto.txt.gz")))
h=f.readline().rstrip("\n").split("\t"); ix={k:i for i,k in enumerate(h)}
out=open(W+"/out/hei/bbj_hei_prscs.txt.tmp","w"); out.write("SNP\tA1\tA2\tBETA\tP\n")
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
out.close(); import os; os.replace(W+"/out/hei/bbj_hei_prscs.txt.tmp", W+"/out/hei/bbj_hei_prscs.txt")
chosen=[]; merged=[]
res=dict(n_records=n, n_sig=len(sig), n_leads=len(chosen), n_loci=len(merged), N_est_median=float(np.median(nest)), N_est_iqr=[float(np.percentile(nest,25)),float(np.percentile(nest,75))],
         loci_kb_median=None)
json.dump(res, open(W+"/out/hei/prscs_prep_summary.json","w"), indent=1); print(json.dumps(res))
