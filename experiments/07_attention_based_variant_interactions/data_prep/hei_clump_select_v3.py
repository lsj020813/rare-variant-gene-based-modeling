#!/usr/bin/env python
import os as _os
if not _os.environ.get("PROJECT_ROOT", "").strip():
    raise ValueError("PROJECT_ROOT must be set and nonblank")
import gzip, zipfile, io, json, os, glob, subprocess, collections, random, csv, numpy as np
def csv_rows(f): return csv.DictReader(open(f),delimiter='\t')
from multiprocessing import Pool
W=(_os.environ["PROJECT_ROOT"] + '/work/prs'); R=(_os.environ["PROJECT_ROOT"] + '/work/ref'); B=_os.environ.get("BCFTOOLS_BIN", "bcftools")
PL=_os.environ.get("PLINK2_BIN", "plink2")
O=W+"/out/hei"; S=O+"/sel"; C=O+"/clump"; os.makedirs(S,exist_ok=True); os.makedirs(C,exist_ok=True); CAP=256
def split_bbj():
    if os.path.exists(C+"/split.done"): return
    Z=R+"/bbj/hum0197.v3.BBJ.Hei.v1.zip"
    f=io.TextIOWrapper(gzip.GzipFile(fileobj=zipfile.ZipFile(Z).open("hum0197.v3.BBJ.Hei.v1/GWASsummary_Height_Japanese_SakaueKanai2020.auto.txt.gz")))
    h=f.readline().rstrip("\n").split("\t"); ix={k:i for i,k in enumerate(h)}
    outs={str(n):open(f"{C}/bbj_chr{n}.tsv","w") for n in range(1,23)}
    for o in outs.values(): o.write("ID\tP\n")
    nsig=0
    for line in f:
        p=line.rstrip("\n").split("\t"); ch=p[ix["CHR"]]; bp=int(p[ix["BP"]])
        if ch not in outs or (ch=="6" and 25_000_000<=bp<=34_000_000): continue
        a1=p[ix["ALLELE1"]]; a0=p[ix["ALLELE0"]]; pv=p[ix["P_BOLT_LMM_INF"]]
        if float(pv)<5e-8: nsig+=1
        outs[ch].write(f"{ch}:{bp}:{a0}:{a1}\t{pv}\n{ch}:{bp}:{a1}:{a0}\t{pv}\n")
    for o in outs.values(): o.close()
    json.dump({"n_sig_nonMHC":nsig},open(C+"/split.done","w"))
def keep_file():
    kf=C+"/ld_keep.txt"
    if os.path.exists(kf): return kf
    sp=[l.rstrip("\n").split("\t") for l in open(W+"/private/hei/split_hei.tsv")][1:]
    trn=[a for a,b in sp if b=="train"]; random.Random(20260927).shuffle(trn)
    open(kf,"w").write("\n".join(trn[:10000])+"\n"); return kf

def susie_cs():
    f=C+"/susie_cs.tsv"
    if os.path.exists(f): return
    Z=R+"/bbj_fm/hum0197.v5.finemap.Hei.v1.zip"
    g=io.TextIOWrapper(gzip.GzipFile(fileobj=zipfile.ZipFile(Z).open("BBJ.Height.Kanai2021.SuSiE.tsv.gz")))
    h=g.readline().rstrip("\n").split("\t"); ix={k:i for i,k in enumerate(h)}
    with open(f+".tmp","w") as o:
        o.write("chr\tpos19\ta1\ta2\tregion\tcs_id\tpip\n")
        for l in g:
            p=l.rstrip("\n").split("\t")
            if p[ix["cs_id"]]=="-1": continue
            o.write("\t".join(p[ix[k]] for k in ("chromosome","position","allele1","allele2","region","cs_id","pip"))+"\n")
    os.replace(f+".tmp",f)
def ld_proxies(args):
    ch,leads=args; o=f"{C}/ld_chr{ch}"
    if not leads: return ch,{}
    if not os.path.exists(o+".vcor"):
        open(o+".leads","w").write("\n".join(k for k,_ in leads)+"\n")
        subprocess.run([PL,"--pfile",f"{W}/geno_hm3/chr{ch}","--keep",C+"/ld_keep.txt","--r2-unphased","--ld-snp-list",o+".leads",
                        "--ld-window-kb","1000","--ld-window","999999","--ld-window-r2","0.5","--threads","2","--memory","6000","--out",o],
                       stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    px=collections.defaultdict(set)
    if os.path.exists(o+".vcor"):
        with open(o+".vcor") as g:
            h=g.readline().lstrip("#").split(); ia=h.index("ID_A"); ib=h.index("ID_B")
            for l in g:
                a=l.split(); px[a[ia]].add(a[ib])
    return ch,px
def clump(ch):
    o=f"{C}/chr{ch}"
    if not os.path.exists(o+".clumps"):
        subprocess.run([PL,"--pfile",f"{W}/geno_hm3/chr{ch}","--keep",C+"/ld_keep.txt","--clump",f"{C}/bbj_chr{ch}.tsv","--clump-p1","5e-8","--clump-p2","1e-2",
                        "--clump-r2","0.1","--clump-kb","1000","--clump-id-field","ID","--clump-p-field","P","--threads","2","--memory","6000","--out",o],
                       stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    leads=[]
    if os.path.exists(o+".clumps"):
        with open(o+".clumps") as g:
            h=g.readline().lstrip("#").split(); ii=h.index("ID"); ip=h.index("P")
            for l in g:
                a=l.split(); leads.append((a[ii],float(a[ip])))
    return ch,leads
def genes_elements(ch):
    genes=[]; cds=[]
    with gzip.open(f"{R}/deductive/gencode.nochr.gtf.gz","rt") as g:
        for l in g:
            if l[0]=="#": continue
            a=l.split("\t",9)
            if a[0].replace("chr","")!=ch: continue
            if a[2]=="CDS": cds.append((int(a[3]),int(a[4])))
            elif a[2]=="gene" and 'gene_type "protein_coding"' in a[8]:
                genes.append((a[8].split('gene_name "')[1].split('"')[0],int(a[3]),int(a[4])))
    el=collections.defaultdict(list)
    for fn in sorted(glob.glob(R+"/re2g/*.bed.gz")):
        with gzip.open(fn,"rt") as g:
            h=g.readline().rstrip("\n").split("\t"); it=h.index("TargetGene"); ia=h.index("ABC.Score")
            for l in g:
                a=l.rstrip("\n").split("\t")
                if a[0].replace("chr","")!=ch: continue
                sc=float(a[ia]) if a[ia] not in ("","NA") else 0.0
                el[a[it]].append((int(a[1]),int(a[2]),sc))
    cds.sort(); return genes,el,cds
def per_chr(args):
    ch,leads,px,csl=args
    if not leads: return ch,[],[],[],{"unc":[],"src":[]}
    ids={k for k,_ in leads}
    for k in list(ids): ids|=px.get(k,set())
    for k,L in csl.items():
        for p19,a1,a2 in L: ids.add(f"{ch}:{p19}:{a1}:{a2}"); ids.add(f"{ch}:{p19}:{a2}:{a1}")
    l38={}
    with gzip.open(f"{R}/lift38_keyed/chr{ch}.keyed38.vcf.gz","rt") as g:
        for l in g:
            if l[0]=="#": continue
            a=l.split("\t",3); kk=a[2][3:] if a[2].startswith("chr") else a[2]
            if kk in ids and a[0].replace("chr","")==ch: l38[kk]=int(a[1])
    genes,el,cds=genes_elements(ch)
    gb={nm:(s,e) for nm,s,e in genes}
    cs=np.array([x[0] for x in cds]); ce=np.maximum.accumulate(np.array([x[1] for x in cds]))
    def in_cds(p):
        i=np.searchsorted(cs,p,side="right")-1; return i>=0 and ce[i]>=p
    link=collections.defaultdict(set); lead_rows=[]; n_nolift=0; SRC={}
    for k,pv in leads:
        if k not in l38: n_nolift+=1; continue
        p=l38[k]; lead_rows.append((ch,k.split(":")[1],pv,k,p))
        lv={p}|{l38[q] for q in px.get(k,()) if q in l38}
        for p19,a1,a2 in csl.get(k,[]):
            for q in (f"{ch}:{p19}:{a1}:{a2}",f"{ch}:{p19}:{a2}:{a1}"):
                if q in l38: lv.add(l38[q])
        SRC[k]=(len(px.get(k,())),len(csl.get(k,[])),len(lv))
        for nm,s,e in genes:
            if any(s-3000<=q<=e+3000 for q in lv): link[k].add(nm)
        for nm,L in el.items():
            if nm in gb and any(s<=q<=e for s,e,_ in L for q in lv): link[k].add(nm)
        if genes:
            nm=min(genes,key=lambda x:0 if x[1]<=p<=x[2] else min(abs(p-x[1]),abs(p-x[2])))[0]; link[k].add(nm)
    sel_genes=sorted(set().union(*link.values())) if link else []
    iv={}
    for nm in sel_genes:
        s,e=gb[nm]; iv[nm]=[(max(1,s-3000),e+3000,0.0)]+[(a,b,sc) for a,b,sc in el.get(nm,[])]
    regs=sorted({(a,b) for L in iv.values() for a,b,_ in L}); rf=f"{S}/q_chr{ch}.regions"
    open(rf,"w").write("".join(f"chr{ch}\t{a}\t{b}\n" for a,b in regs))
    out=subprocess.run([B,"query","-R",rf,"-i","INFO/R2>=0.3 && INFO/MAF>=0.0008","-f","%CHROM\t%POS\t%ID\n",f"{R}/lift38_keyed/chr{ch}.keyed38.vcf.gz"],capture_output=True,text=True).stdout
    os.remove(rf)
    V={}
    for l in out.splitlines():
        c_,p,k=l.split("\t"); k=k[3:] if k.startswith("chr") else k
        if c_.replace("chr","")!=ch or not k.startswith(ch+":"): continue
        p=int(p)
        if not in_cds(p): V[k]=p
    items=sorted(V.items(),key=lambda x:x[1]); pos=np.array([p for _,p in items]); keys=[k for k,_ in items]
    grows=[]; gtok={}; UNC={}
    for nm in sel_genes:
        s,e=gb[nm]; d={}
        for a,b,sc in iv[nm]:
            for i in range(np.searchsorted(pos,a),np.searchsorted(pos,b,side="right")):
                if sc>d.get(i,-1): d[i]=sc
        dist=lambda i: 0 if s<=pos[i]<=e else min(abs(pos[i]-s),abs(pos[i]-e))
        full=sorted(d,key=lambda i:(-d[i],dist(i),pos[i])); UNC[nm]=len(full); order=full[:CAP]; gtok[nm]={keys[i] for i in order}
        for r_,i in enumerate(order): grows.append((nm,ch,keys[i],int(pos[i]),round(d[i],6),r_))
    drows=[]
    for (c_,bp19,pv,k,p38) in lead_rows:
        n=len(set().union(*[gtok[g] for g in link[k]])) if link[k] else 0
        if n==0: continue
        o2=subprocess.run([B,"query","-r",f"chr{ch}:{max(1,p38-1000000)}-{p38+1000000}","-i","INFO/R2>=0.3 && INFO/MAF>=0.0008","-f","%CHROM\t%POS\t%ID\n",f"{R}/lift38_keyed/chr{ch}.keyed38.vcf.gz"],capture_output=True,text=True).stdout
        cand=[]
        for l in o2.splitlines():
            c2,pp,kk=l.split("\t"); kk=kk[3:] if kk.startswith("chr") else kk
            if c2.replace("chr","")==ch and kk.startswith(ch+":") and not in_cds(int(pp)): cand.append((abs(int(pp)-p38),int(pp),kk))
        cand.sort()
        for _,pp,kk in cand[:n]: drows.append((k,kk,pp))
    allk=sorted({r_[2] for r_ in grows}|{r_[1] for r_ in drows},key=lambda x:int(x.split(":")[1]))
    open(f"{S}/extract_chr{ch}.keys","w").write("\n".join(allk)+("\n" if allk else ""))
    return ch,lead_rows,grows,drows,dict(n_leads=len(leads),n_nolift=n_nolift,n_genes=len(sel_genes),links=sum(len(v) for v in link.values()),unc=list(UNC.values()),src=list(SRC.values()))
if __name__=="__main__":
    split_bbj(); keep_file()
    chs=[str(n) for n in range(1,23)]
    with Pool(6) as p: cl=p.map(clump,chs)
    susie_cs()
    with Pool(6) as p: pxl=dict(p.map(ld_proxies,cl))
    CS=collections.defaultdict(list)
    for row in csv_rows(C+"/susie_cs.tsv"): CS[(row["chr"],row["region"],row["cs_id"])].append((int(row["pos19"]),row["a1"],row["a2"],float(row["pip"])))
    leadpos={ch:[(int(k.split(":")[1]),k) for k,_ in L] for ch,L in cl}
    asg=collections.defaultdict(lambda: collections.defaultdict(list)); n_cs=0; n_cs_asg=0
    for (ch,reg,cid),vs in CS.items():
        n_cs+=1; top=max(vs,key=lambda x:x[3])[0]; lp=leadpos.get(ch,[])
        if not lp: continue
        d,k=min((abs(top-p),k) for p,k in lp)
        if d<=500000: n_cs_asg+=1; asg[ch][k]+= [(p,a1,a2) for p,a1,a2,_ in vs]
    CSSUM=dict(n_cs=n_cs,n_cs_assigned=n_cs_asg)
    with Pool(6) as p: res=p.map(per_chr,[(ch,L,pxl.get(ch,{}),dict(asg[ch])) for ch,L in cl])
    with open(O+"/lead_snps.tsv","w") as g:
        for ch,lr,_,_,_ in res:
            for c_,bp19,pv,k,p38 in lr: g.write(f"{c_}\t{bp19}\t{pv}\t{k}\t{p38}\n")
    with open(S+"/gene_tokens.tsv.tmp","w") as g:
        g.write("gene\tchr\tkey19\tpos38\tabc\trank\n")
        for _,_,gr,_,_ in res:
            for x in gr: g.write("\t".join(map(str,x))+"\n")
    with open(S+"/dist_tokens.tsv.tmp","w") as g:
        g.write("lead\tkey19\tpos38\n")
        for _,_,_,dr,_ in res:
            for x in dr: g.write("\t".join(map(str,x))+"\n")
    os.replace(S+"/gene_tokens.tsv.tmp",S+"/gene_tokens.tsv"); os.replace(S+"/dist_tokens.tsv.tmp",S+"/dist_tokens.tsv")
    gr=[x for r_ in res for x in r_[2]]; dr=[x for r_ in res for x in r_[3]]
    ntok=collections.Counter(x[0] for x in gr); nt=np.array(list(ntok.values())) if ntok else np.array([0])
    uk={x[2] for x in gr}; ud={x[1] for x in dr}
    srcs=[x for r_ in res for x in r_[4].get('src',[])]
    summ=dict(rule="C+LDproxy0.5+SuSiE_CS",cs=CSSUM,leads_with_proxy=sum(1 for a in srcs if a[0]>0),leads_with_cs=sum(1 for a in srcs if a[1]>0),link_variants_per_lead_median=float(np.median([a[2] for a in srcs])) if srcs else 0,n_sig_nonMHC=json.load(open(C+"/split.done"))["n_sig_nonMHC"],indep_leads=sum(r_[4].get("n_leads",0) for r_ in res),
              leads_not_lifted=sum(r_[4].get("n_nolift",0) for r_ in res),genes=len(ntok),gene_links=sum(r_[4].get("links",0) for r_ in res),
              tokens_with_dup=len(gr),unique_ann=len(uk),unique_dist=len(ud),overlap=len(uk&ud),unique_to_extract=len(uk|ud),
              tokens_per_gene_q=[float(np.percentile(nt,q)) for q in (0,25,50,75,100)],genes_at_cap=int((nt>=CAP).sum()),
              uncapped_tokens_per_gene_q=[float(np.percentile(np.array([u for r_ in res for u in r_[4].get('unc',[])] or [0]),q)) for q in (0,25,50,75,90,100)],frac_tokens_abc_gt0=round(sum(x[4]>0 for x in gr)/max(1,len(gr)),4))
    json.dump(summ,open(S+"/selection_summary.json","w"),indent=1); print(json.dumps(summ))
