#!/usr/bin/env python
import os as _os
if not _os.environ.get("PROJECT_ROOT", "").strip():
    raise ValueError("PROJECT_ROOT must be set and nonblank")
import gzip, json, os, glob, subprocess, collections, numpy as np
from multiprocessing import Pool
W=(_os.environ["PROJECT_ROOT"] + '/work/prs'); R=(_os.environ["PROJECT_ROOT"] + '/work/ref'); B=_os.environ.get("BCFTOOLS_BIN", "bcftools")
O=W+"/out/hei/sel"; os.makedirs(O,exist_ok=True); CAP=64
leads=collections.defaultdict(list)
for l in open(W+"/out/hei/lead_snps.tsv"):
    ch,bp,_=l.split(); leads[ch].append(int(bp))
def per_chr(ch):
    loci=[l.split() for l in open(W+"/out/hei/loci.bed") if l.split()[0]==ch]
    var=[]
    for _,s,e,L in loci:
        out=subprocess.run([B,"query","-r",f"{ch}:{s}-{e}","-i","INFO/R2>=0.3 && INFO/MAF>=0.0008","-f","%POS\t%REF\t%ALT\n",f"{R}/orig_index/chr{ch}.vcf.gz"],capture_output=True,text=True).stdout
        for x in out.splitlines():
            a=x.split("\t"); var.append((L,f"{ch}:{a[0]}:{a[1]}:{a[2]}",int(a[0])))
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
                genes.append((a[8].split('gene_name "')[1].split('"')[0],int(a[3]),int(a[4])))
    cds.sort(); cs=np.array([x[0] for x in cds]); ce=np.maximum.accumulate(np.array([x[1] for x in cds]))
    def in_cds(p):
        i=np.searchsorted(cs,p,side="right")-1; return i>=0 and ce[i]>=p
    V=[(L,k,p19,k2p[k]) for L,k,p19 in var if k in k2p and not in_cds(k2p[k])]
    V.sort(key=lambda x:x[3]); pos=np.array([v[3] for v in V])
    gb={nm:(s,e) for nm,s,e in genes}
    score=collections.defaultdict(dict)
    for nm,s,e in genes:
        for i in range(np.searchsorted(pos,s-3000),np.searchsorted(pos,e+3000,side="right")): score[nm].setdefault(i,0.0)
    for fn in sorted(glob.glob(R+"/re2g/*.bed.gz")):
        with gzip.open(fn,"rt") as g:
            h=g.readline().rstrip("\n").split("\t"); it=h.index("TargetGene"); ia=h.index("ABC.Score")
            for l in g:
                a=l.rstrip("\n").split("\t")
                if a[0].replace("chr","")!=ch: continue
                s,e=int(a[1]),int(a[2]); sc=float(a[ia]) if a[ia] not in ("","NA") else 0.0
                d=score[a[it]]
                for i in range(np.searchsorted(pos,s),np.searchsorted(pos,e,side="right")):
                    if sc>d.get(i,-1): d[i]=sc
    rows=[]
    for g_,d in score.items():
        if not d: continue
        s,e=gb.get(g_,(None,None))
        def dist(i):
            if s is None: return 10**9
            p=V[i][3]; return 0 if s<=p<=e else min(abs(p-s),abs(p-e))
        order=sorted(d, key=lambda i:(-d[i],dist(i),V[i][3]))[:CAP]
        Lc=collections.Counter(V[i][0] for i in order).most_common(1)[0][0]
        for r,i in enumerate(order): rows.append((g_,Lc,V[i][1],V[i][3],round(d[i],6),r))
    ann_by_locus=collections.defaultdict(set)
    for g_,Lc,k,p38,sc,r in rows: ann_by_locus[Lc].add(k)
    drows=[]
    for _,s,e,L in loci:
        n=len(ann_by_locus.get(L,()))
        if n==0: continue
        ld=[b for b in leads[ch] if int(s)<=b<=int(e)] or [(int(s)+int(e))//2]
        cand=sorted({(k,p19) for LL,k,p19,_ in V if LL==L}, key=lambda kp:(min(abs(kp[1]-b) for b in ld),kp[1]))[:n]
        for k,p19 in cand: drows.append((L,k,p19))
    keys=sorted({r[2] for r in rows}|{r[1] for r in drows}, key=lambda k:int(k.split(":")[1]))
    with open(f"{O}/extract_chr{ch}.keys","w") as g: g.write("\n".join(keys)+"\n")
    return ch,rows,drows,len(V)
if __name__=="__main__":
    chs=sorted({l.split()[0] for l in open(W+"/out/hei/loci.bed")},key=int)
    with Pool(8) as p: res=p.map(per_chr,chs)
    with open(O+"/gene_tokens.tsv.tmp","w") as g:
        g.write("gene\tlocus\tkey19\tpos38\tabc\trank\n")
        for _,rows,_,_ in res:
            for x in rows: g.write("\t".join(map(str,x))+"\n")
    with open(O+"/dist_tokens.tsv.tmp","w") as g:
        g.write("locus\tkey19\tpos19\n")
        for _,_,dr,_ in res:
            for x in dr: g.write("\t".join(map(str,x))+"\n")
    os.replace(O+"/gene_tokens.tsv.tmp",O+"/gene_tokens.tsv"); os.replace(O+"/dist_tokens.tsv.tmp",O+"/dist_tokens.tsv")
    allrows=[x for _,rows,_,_ in res for x in rows]; dr=[x for _,_,d,_ in res for x in d]
    uk=len({x[2] for x in allrows}); ud=len({x[1] for x in dr}); uall=len({x[2] for x in allrows}|{x[1] for x in dr})
    summ=dict(cap=CAP,genes=len({x[0] for x in allrows}),tokens_with_dup=len(allrows),unique_ann=uk,unique_dist=ud,unique_union_to_extract=uall,
              overlap_ann_dist=uk+ud-uall,frac_tokens_with_abc_gt0=round(sum(x[4]>0 for x in allrows)/len(allrows),4))
    json.dump(summ,open(O+"/selection_summary.json","w"),indent=1); print(json.dumps(summ))
