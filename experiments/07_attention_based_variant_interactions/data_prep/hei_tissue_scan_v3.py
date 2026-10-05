#!/usr/bin/env python
import os as _os
if not _os.environ.get("PROJECT_ROOT", "").strip():
    raise ValueError("PROJECT_ROOT must be set and nonblank")
import os, glob, gzip, json, subprocess, numpy as np, pyarrow.parquet as pq
from multiprocessing import Pool
W=(_os.environ["PROJECT_ROOT"] + '/work'); R=W+"/ref"; O=W+"/prs/out/hei/tissue_scan_v3"; os.makedirs(O,exist_ok=True)
B=_os.environ.get("BCFTOOLS_BIN", "bcftools"); rng=np.random.default_rng(20260927); NPERM=1000
MAFB=np.array([0.0008,0.005,0.01,0.05,0.1,0.2,0.51]); DB=np.array([0,1e3,5e3,2e4,5e4,1e5,2.5e5,1e6,1e12])
tss={}; cds={}
with gzip.open(R+"/deductive/gencode.nochr.gtf.gz","rt") as g:
    for l in g:
        if l[0]=="#": continue
        p=l.split("\t",9)
        if 'gene_type "protein_coding"' not in p[8]: continue
        ch=p[0]
        if p[2]=="gene": tss.setdefault(ch,[]).append(int(p[3]) if p[6]=="+" else int(p[4]))
        elif p[2]=="CDS": cds.setdefault(ch,[]).append((int(p[3]),int(p[4])))
def merge(iv):
    iv=sorted(iv); s=[];e=[]
    for a,b in iv:
        if s and a<=e[-1]: e[-1]=max(e[-1],b)
        else: s.append(a); e.append(b)
    return np.array(s,np.int64),np.array(e,np.int64)
tss={k:np.sort(np.array(v,np.int64)) for k,v in tss.items()}; cdsm={k:merge(v) for k,v in cds.items()}
def incds(ch,p):
    if ch not in cdsm: return np.zeros(len(p),bool)
    s,e=cdsm[ch]; i=np.searchsorted(s,p,"right")-1; return (i>=0)&(p<=e[np.maximum(i,0)])
def tssdist(ch,p):
    t=tss.get(ch); 
    if t is None: return np.full(len(p),1e12)
    i=np.searchsorted(t,p); a=np.abs(p-t[np.clip(i-1,0,len(t)-1)]); b=np.abs(t[np.clip(i,0,len(t)-1)]-p); return np.minimum(a,b)
cs=[l.rstrip("\n").split("\t") for l in open(W+"/prs/out/hei/clump/susie_cs.tsv")][1:]
want={}
for c,p19,a1,a2,reg,cid,pip in cs:
    for k in (f"chr{c}:{p19}:{a1}:{a2}",f"chr{c}:{p19}:{a2}:{a1}"): want[k]=(reg+"|"+cid,float(pip))
def uni(N):
    q=subprocess.run([B,"query","-i","INFO/R2>=0.3 && INFO/MAF>=0.0008","-f","%CHROM\t%POS\t%ID\t%INFO/MAF\n",f"{R}/lift38_keyed/chr{N}.keyed38.vcf.gz"],capture_output=True,text=True,check=True).stdout
    P=[];M=[];hit=[]
    for l in q.splitlines():
        c38,p38,k,m=l.split("\t")
        if c38!=f"chr{N}": continue
        P.append(int(p38)); M.append(float(m))
        if k in want: hit.append((want[k][0],want[k][1],int(p38),float(m)))
    P=np.array(P,np.int64); M=np.array(M); ch=str(N); nc=~incds(ch,P)
    return N,P[nc],M[nc],tssdist(ch,P[nc]),[h for h in hit if not incds(ch,np.array([h[2]]))[0]]
with Pool(8) as pl: U=pl.map(uni,range(1,23))
pool_ch=[];pool_p=[];pool_s=[]; csv={}
for N,P,M,D,hits in U:
    s=(np.digitize(M,MAFB)-1)*10+(np.digitize(D,DB)-1); pool_ch.append(np.full(len(P),N,np.int16)); pool_p.append(P); pool_s.append(s)
    for key,pip,p38,m in hits: csv.setdefault(key,[]).append((N,p38,pip,m))
pool_ch=np.concatenate(pool_ch); pool_p=np.concatenate(pool_p); pool_s=np.concatenate(pool_s)
bys={s:np.where(pool_s==s)[0] for s in np.unique(pool_s)}
CS=[]
for key,v in csv.items():
    v=sorted(v,key=lambda x:-x[2]); N,lp,_,lm=v[0]
    ls=(np.digitize([lm],MAFB)[0]-1)*10+(np.digitize(tssdist(str(N),np.array([lp])),DB)[0]-1)
    CS.append(dict(key=key,chr=N,pos=np.array([x[1] for x in v]),pip=np.array([x[2] for x in v]),off=np.array([x[1]-lp for x in v]),stratum=int(ls)))
CS=[c for c in CS if c["stratum"] in bys]
info=dict(n_cs_raw=len(set(x[0] for x in want.values())),n_cs_used=len(CS),n_var_used=int(sum(len(c["pos"]) for c in CS)),pool=int(len(pool_p)))
print(json.dumps(info),flush=True)
anc=np.stack([rng.choice(bys[c["stratum"]],NPERM) for c in CS],1)
T={}
for f in sorted(glob.glob(R+"/re2g_all/ot_e2g/*.parquet")):
    t=pq.read_table(f,columns=["biosampleName","chromosome","start","end"]).to_pandas()
    for (b,ch),g in t.groupby(["biosampleName","chromosome"]):
        T.setdefault(b,{}).setdefault(ch,[]).append(np.stack([g.start.values,g.end.values],1))
names=sorted(T); print("tissues",len(names),flush=True)
def inel(iv,ch,p):
    if ch not in iv: return np.zeros(len(p),bool)
    s,e=iv[ch]; i=np.searchsorted(s,p,"left")-1; return (i>=0)&(p<=e[np.maximum(i,0)])
def merge_arr(x):
    x=x[np.argsort(x[:,0])]; s_=[];e_=[]
    for u,v in x:
        if s_ and u<=e_[-1]:
            if v>e_[-1]: e_[-1]=v
        else: s_.append(u); e_.append(v)
    return np.array(s_,np.int64),np.array(e_,np.int64)
def one(name):
    iv={ch:merge_arr(np.concatenate(v)) for ch,v in T[name].items()}
    obs=np.mean([(c["pip"]*inel(iv,str(c["chr"]),c["pos"])).sum()/c["pip"].sum() for c in CS])
    null=np.zeros(NPERM)
    for j,c in enumerate(CS):
        a_=anc[:,j]; ch=pool_ch[a_]; pos=pool_p[a_][:,None]+c["off"][None,:]; w=c["pip"]/c["pip"].sum(); m=np.zeros(pos.shape,bool)
        for N in np.unique(ch):
            r_=ch==N; m[r_]=inel(iv,str(N),pos[r_].ravel()).reshape(r_.sum(),-1)
        null+=(m*w[None,:]).sum(1)
    null/=len(CS)
    return name,float(obs),null,int(sum(len(v[0]) for v in iv.values()))
with Pool(8) as pl: out=pl.map(one,names,chunksize=4)
EPS=1e-6
L=np.log(np.array([o[1] for o in out])+EPS); LN=np.log(np.stack([o[2] for o in out])+EPS)
D=L-L.mean(); DN=LN-LN.mean(0,keepdims=True)
p=(1+(DN>=D[:,None]).sum(1))/(NPERM+1)
res=[dict(tissue=o[0],S=o[1],null_mean=float(o[2].mean()),null_sd=float(o[2].std()),enrich=float(o[1]/max(o[2].mean(),1e-12)),
          D=float(D[i]),rel_fold=float(np.exp(D[i]-DN[i].mean())),p=float(p[i]),n_elem=o[3]) for i,o in enumerate(out)]
o_=np.argsort(p); q=np.empty(len(p)); m_=len(p); q[o_]=np.minimum.accumulate((p[o_]*m_/np.arange(1,m_+1))[::-1])[::-1]
for r,qq in zip(res,q): r["fdr"]=float(min(qq,1.0)); r["seed"]=bool(qq<0.05 and r["D"]>0)
res.sort(key=lambda r:(r["p"],-r["D"]))
with open(O+"/tissue_scan.tsv.tmp","w") as g:
    g.write("tissue\tS\tnull_mean\tnull_sd\tenrich\tD\trel_fold\tp\tfdr\tseed\tn_elem\n")
    for r in res: g.write("\t".join(str(r[k]) for k in ("tissue","S","null_mean","null_sd","enrich","D","rel_fold","p","fdr","seed","n_elem"))+"\n")
os.replace(O+"/tissue_scan.tsv.tmp",O+"/tissue_scan.tsv")
seeds=[r["tissue"] for r in res if r["seed"]]; info.update(n_seeds=len(seeds),fallback_to_literature=len(seeds)<3,seeds=seeds)
json.dump(info,open(O+"/tissue_scan_summary.json","w"),indent=1); print(json.dumps(info)[:3000],flush=True)
