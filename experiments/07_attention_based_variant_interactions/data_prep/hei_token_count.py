#!/usr/bin/env python
import os as _os
if not _os.environ.get("PROJECT_ROOT", "").strip():
    raise ValueError("PROJECT_ROOT must be set and nonblank")
import gzip, zipfile, io, json, os, subprocess, sys, bisect, collections, numpy as np
from multiprocessing import Pool
W=(_os.environ["PROJECT_ROOT"] + '/work/prs'); R=(_os.environ["PROJECT_ROOT"] + '/work/ref'); B=_os.environ.get("BCFTOOLS_BIN", "bcftools")
O=W+"/out/hei"; os.makedirs(O,exist_ok=True)
def make_loci():
    if os.path.exists(O+"/loci.bed"): return
    Z=R+"/bbj/hum0197.v3.BBJ.Hei.v1.zip"
    f=io.TextIOWrapper(gzip.GzipFile(fileobj=zipfile.ZipFile(Z).open("hum0197.v3.BBJ.Hei.v1/GWASsummary_Height_Japanese_SakaueKanai2020.auto.txt.gz")))
    h=f.readline().rstrip("\n").split("\t"); ix={k:i for i,k in enumerate(h)}; sig=[]
    for line in f:
        p=line.split("\t"); pv=float(p[ix["P_BOLT_LMM_INF"]]); ch=p[ix["CHR"]]; bp=int(p[ix["BP"]])
        if pv<5e-8 and not (ch=="6" and 25_000_000<=bp<=34_000_000): sig.append((pv,ch,bp))
    sig.sort(); chosen=[]
    for pv,ch,bp in sig:
        if all(not (c==ch and abs(b-bp)<=250000) for _,c,b in chosen): chosen.append((pv,ch,bp))
    wins=sorted(((ch,max(0,bp-250000),bp+250000) for _,ch,bp in chosen), key=lambda x:(int(x[0]),x[1])); m=[]
    for ch,s,e in wins:
        if m and m[-1][0]==ch and s<=m[-1][2]: m[-1][2]=max(m[-1][2],e)
        else: m.append([ch,s,e])
    with open(O+"/loci.bed.tmp","w") as g:
        for i,(ch,s,e) in enumerate(m): g.write(f"{ch}\t{s}\t{e}\tH{i:03d}\n")
    with open(O+"/lead_snps.tsv","w") as g:
        for pv,ch,bp in sorted(chosen,key=lambda x:(int(x[1]),x[2])): g.write(f"{ch}\t{bp}\t{pv}\n")
    os.replace(O+"/loci.bed.tmp",O+"/loci.bed")
def per_chr(ch):
    loci=[l.split() for l in open(O+"/loci.bed") if l.split()[0]==ch]
    if not loci: return None
    var=[]
    for _,s,e,L in loci:
        out=subprocess.run([B,"query","-r",f"{ch}:{s}-{e}","-i","INFO/R2>=0.3 && INFO/MAF>=0.0008","-f","%POS\t%REF\t%ALT\t%INFO/MAF\n",f"{R}/orig_index/chr{ch}.vcf.gz"],capture_output=True,text=True).stdout
        for l in out.splitlines():
            a=l.split("\t"); var.append((L,f"{ch}:{a[0]}:{a[1]}:{a[2]}",float(a[3])))
    want={k for _,k,_ in var}; k2p={}
    with gzip.open(f"{R}/lift38_keyed/chr{ch}.keyed38.vcf.gz","rt") as g:
        for l in g:
            if l[0]=="#": continue
            a=l.split("\t",3); kk=a[2][3:] if a[2].startswith("chr") else a[2]
            if kk in want and a[0].replace("chr","")==ch: k2p[kk]=int(a[1])
    cds=[]; genes=[]
    with gzip.open(f"{R}/deductive/gencode.nochr.gtf.gz","rt") as g:
        for l in g:
            if l[0]=="#": continue
            a=l.split("\t",9)
            if a[0].replace("chr","")!=ch: continue
            if a[2]=="CDS": cds.append((int(a[3]),int(a[4])))
            elif a[2]=="gene" and 'gene_type "protein_coding"' in a[8]:
                nm=a[8].split('gene_name "')[1].split('"')[0]; genes.append((nm,int(a[3]),int(a[4])))
    cds.sort(); cs=np.array([x[0] for x in cds]); ce=np.maximum.accumulate(np.array([x[1] for x in cds])) if cds else np.array([])
    def in_cds(p):
        if not cds: return False
        i=np.searchsorted(cs,p,side="right")-1
        return i>=0 and ce[i]>=p
    V=[(L,k,maf,k2p[k]) for L,k,maf in var if k in k2p]
    n_lift=len(V); V=[v for v in V if not in_cds(v[3])]; n_nc=len(V)
    V.sort(key=lambda x:x[3]); pos=np.array([v[3] for v in V])
    assign=collections.defaultdict(set); src=collections.Counter(); bios=set()
    for nm,s,e in genes:
        a=np.searchsorted(pos,s-3000); b=np.searchsorted(pos,e+3000,side="right")
        for i in range(a,b): assign[nm].add(i); src["body3kb"]+=1
    import glob
    for fn in glob.glob(R+"/re2g/*.bed.gz"):
        with gzip.open(fn,"rt") as g:
            hdr=g.readline().rstrip("\n").split("\t"); it=hdr.index("TargetGene")
            for l in g:
                a=l.rstrip("\n").split("\t")
                if a[0].replace("chr","")!=ch: continue
                if len(a)>27: bios.add(a[27])
                s,e=int(a[1]),int(a[2]); lo=np.searchsorted(pos,s); hi=np.searchsorted(pos,e,side="right")
                for i in range(lo,hi): assign[a[it]].add(i); src["re2g"]+=1
    anyv=set().union(*assign.values()) if assign else set()
    mult=collections.Counter(); 
    for g_,ss in assign.items():
        for i in ss: mult[i]+=1
    per_gene=[]
    for g_,ss in assign.items():
        if not ss: continue
        Ls=collections.Counter(V[i][0] for i in ss); mafs=[min(V[i][2],1-V[i][2]) for i in ss]
        per_gene.append((g_,ch,Ls.most_common(1)[0][0],len(ss),sum(m>=0.05 for m in mafs),sum(0.01<=m<0.05 for m in mafs),sum(m<0.01 for m in mafs)))
    per_locus=[]
    for _,s,e,L in loci:
        idx=[i for i,v in enumerate(V) if v[0]==L]; asg=[i for i in idx if i in anyv]
        gs={g_ for g_,ss in assign.items() if any(V[i][0]==L for i in ss)}
        per_locus.append((L,ch,int(e)-int(s),sum(1 for v in var if v[0]==L),len(idx),len(asg),len(gs)))
    return dict(ch=ch,n_query=len(var),n_lift=n_lift,n_noncoding=n_nc,n_assigned_unique=len(anyv),n_multi_gene=sum(1 for i,c in mult.items() if c>1),
                per_gene=per_gene,per_locus=per_locus,bios=sorted(bios),src=dict(src))
if __name__=="__main__":
    make_loci()
    chs=sorted({l.split()[0] for l in open(O+"/loci.bed")},key=int)
    with Pool(8) as p: res=[r for r in p.map(per_chr,chs) if r]
    pg=[x for r in res for x in r["per_gene"]]; pl=[x for r in res for x in r["per_locus"]]
    with open(O+"/token_per_gene.tsv","w") as g:
        g.write("gene\tchr\tlocus\tn_tokens\tn_common\tn_low\tn_rare\n"); [g.write("\t".join(map(str,x))+"\n") for x in pg]
    with open(O+"/token_per_locus.tsv","w") as g:
        g.write("locus\tchr\tbp\tn_window\tn_noncoding_lifted\tn_assigned\tn_genes\n"); [g.write("\t".join(map(str,x))+"\n") for x in pl]
    nt=np.array([x[3] for x in pg]); 
    summ=dict(n_loci=len(pl),n_genes=len(pg),window_variants=int(sum(r["n_query"] for r in res)),noncoding_lifted=int(sum(r["n_noncoding"] for r in res)),
              assigned_unique=int(sum(r["n_assigned_unique"] for r in res)),multi_gene_variants=int(sum(r["n_multi_gene"] for r in res)),
              tokens_total_with_dup=int(nt.sum()),tokens_per_gene_q=[float(np.percentile(nt,q)) for q in (0,10,25,50,75,90,99,100)],
              genes_gt256=int((nt>256).sum()),genes_gt512=int((nt>512).sum()),genes_gt1024=int((nt>1024).sum()),
              loci_with_0_genes=int(sum(1 for x in pl if x[6]==0)),biosamples=sorted(set().union(*[set(r["bios"]) for r in res])),
              maf_mix_total={"common":int(sum(x[4] for x in pg)),"low":int(sum(x[5] for x in pg)),"rare":int(sum(x[6] for x in pg))})
    json.dump(summ,open(O+"/token_count_summary.json","w"),indent=1); print(json.dumps(summ))
