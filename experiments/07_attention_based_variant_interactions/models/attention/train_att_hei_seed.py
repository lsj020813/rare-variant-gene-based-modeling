#!/usr/bin/env python
import os as _os
if not _os.environ.get("PROJECT_ROOT", "").strip():
    raise ValueError("PROJECT_ROOT must be set and nonblank")
import sys, os, json, time, csv, gzip, glob, collections, numpy as np, torch, torch.nn as nn, torch.utils.checkpoint
W=(_os.environ["PROJECT_ROOT"] + '/work/prs'); P0=W+"/private/hei"; P=W+"/private/hei_seed"; O=W+"/out/hei_seed"; SEL=W+"/out/hei/sel_seed"
ARM=sys.argv[1]; MODE=sys.argv[2] if len(sys.argv)>2 else "train"
torch.manual_seed(20260927); np.random.seed(20260927)
d=np.load(P0+"/base_hei.npz",allow_pickle=True); r=d["r"].astype(np.float64); prs=d["prs"].astype(np.float64); sp=d["split"]
tr=sp=="train"; va=sp=="valid"; te=sp=="test"; n=len(r)
A=np.column_stack([np.ones(n),prs]); ba=np.linalg.lstsq(A[tr],r[tr],rcond=None)[0]; base_pred=A@ba; rr=(r-base_pred).astype(np.float32)
def r2(y,p): return float(1-((y-p)**2).sum()/((y-y.mean())**2).sum())
CAP=128; M=CAP
toks=collections.defaultdict(list)
for row in csv.DictReader(open(SEL+"/gene_tokens.tsv"),delimiter="\t"):
    if int(row["rank"])<CAP: toks[row["gene"]].append((int(row["rank"]),row["key19"],float(row["abc"])))
genes=sorted(toks); G=len(genes)
U=sorted({k for g in genes for _,k,_ in toks[g]}, key=lambda k:(int(k.split(":")[0]),int(k.split(":")[1]),k)); ui={k:i for i,k in enumerate(U)}
idx=np.zeros((G,M),np.int64); mask=np.ones((G,M),bool)
for gi,g in enumerate(genes):
    for t,(rk,k,_) in enumerate(sorted(toks[g])[:M]): idx[gi,t]=ui[k]; mask[gi,t]=False
meta={}
for f in glob.glob(O+"/meta/meta_chr*.tsv"):
    for row in csv.DictReader(open(f),delimiter="\t"): meta[row["key"]]=row
assert all(k in meta for k in U), "meta missing"
if ARM=="P5h":
    loc={}
    for f in glob.glob(O+"/gt/k_chr*.tsv"):
        N=os.path.basename(f)[5:-4]
        for j,row in enumerate(csv.DictReader(open(f),delimiter="\t")): loc[row["key"]]=(N,j)
    assert all(k in loc for k in U), f"gt missing {sum(k not in loc for k in U)}"
    H=np.zeros((n,len(U),2),np.int8); byc=collections.defaultdict(list)
    for k in U: byc[loc[k][0]].append((loc[k][1],ui[k]))
    for N,lst in byc.items():
        Hc=np.load(f"{P}/gt/h_chr{N}.npy",mmap_mode="r"); cols=np.array([a for a,_ in lst]); dst=np.array([b for _,b in lst])
        for s in range(0,len(cols),4096): H[:,dst[s:s+4096],:]=Hc[:,cols[s:s+4096],:]
        del Hc
    pf=H[tr].mean(axis=(0,2)).astype(np.float32); sf=np.sqrt(np.maximum(pf*(1-pf),1e-4)).astype(np.float32)
    Pt=torch.tensor(pf); Sf=torch.tensor(sf)
else:
    loc={}
    for f in glob.glob(O+"/meta/meta_chr*.tsv"):
        N=os.path.basename(f)[8:-4]
        for j,row in enumerate(csv.DictReader(open(f),delimiter="\t")): loc[row["key"]]=(N,j)
    D=np.zeros((n,len(U)),np.float16); byc=collections.defaultdict(list)
    for k in U: byc[loc[k][0]].append((loc[k][1],ui[k]))
    for N,lst in byc.items():
        Gm=np.load(f"{P}/geno_chr{N}.npy",mmap_mode="r"); cols=np.array([a for a,_ in lst]); dst=np.array([b for _,b in lst])
        for s in range(0,len(cols),4096):
            X=np.asarray(Gm[:,cols[s:s+4096]],dtype=np.float32); mu=X[tr].mean(0); sd=X[tr].std(0); sd[sd<1e-6]=1; D[:,dst[s:s+4096]]=((X-mu)/sd).astype(np.float16)
        del Gm
F={}; cols_f=None
for f in sorted(glob.glob(O+"/feat/feat_chr*.tsv.gz")):
    with gzip.open(f,"rt") as g:
        h=g.readline().rstrip("\n").split("\t")
        if cols_f is None: cols_f=[c for c in h if c not in ("key","chr","pos38","domain") and "re2g" not in c.lower() and not c.endswith("_target")]
        ix=[h.index(c) for c in cols_f]
        for l in g:
            p_=l.rstrip("\n").split("\t"); F[p_[0]]=[p_[i] for i in ix]
fmiss=sum(k not in F for k in U); assert fmiss==0, f"features missing for {fmiss} tokens"
num=[]; cat=[]
for j,c in enumerate(cols_f):
    vals=[F[k][j] for k in U[:5000] if k in F]
    try: [float(v) for v in vals if v!=""]; num.append(j)
    except ValueError: cat.append(j)
blocks=[]
for j in num:
    v=np.array([float(F[k][j]) if (k in F and F[k][j]!="") else np.nan for k in U],np.float32); m_=np.isnan(v); v[m_]=0; blocks.append(v)
    if m_.any() and not m_.all(): blocks.append(m_.astype(np.float32))
for j in cat:
    for x in sorted({F[k][j] for k in U if k in F})[:30]: blocks.append(np.array([1.0 if (k in F and F[k][j]==x) else 0.0 for k in U],np.float32))
abc={k:0.0 for k in U}
for g in genes:
    for _,k,a in toks[g]: abc[k]=max(abc[k],a)
blocks+=[np.log10(np.maximum(np.array([float(meta[k]["maf_train"]) for k in U]),1e-4)).astype(np.float32),np.array([float(meta[k]["r2"]) for k in U],np.float32),
         np.array([1.0 if meta[k]["typed"] not in ("0","",".") else 0.0 for k in U],np.float32),np.log1p(100*np.array([abc[k] for k in U],np.float32))]
Fk=np.column_stack(blocks).astype(np.float32); mu=Fk.mean(0); sd=Fk.std(0); sd[sd<1e-6]=1; Fk=(Fk-mu)/sd
S=Fk[idx]; S[mask]=0; S=np.concatenate([S,np.tile(np.arange(M)/M,(G,1)).astype(np.float32)[...,None]],-1); fs=S.shape[-1]
dev=torch.device("cuda"); St=torch.tensor(S,device=dev); Mt=torch.tensor(mask,device=dev); IDX=idx.reshape(-1)
NS=2 if ARM=="P5h" else 1; SCALE=float(np.sqrt(G*M*NS))
print(json.dumps(dict(arm=ARM,mode=MODE,genes=G,M=M,U=len(U),fs=fs,n=n,strands=NS,scale=SCALE)),flush=True)
class Net(nn.Module):
    def __init__(s,d=32,h=4,drop=0.1,mix=True):
        super().__init__()
        s.fe=nn.Sequential(nn.Linear(fs,d),nn.GELU(),nn.Linear(d,d)); s.fd=nn.Linear(fs,d); s.le=nn.Embedding(G,d); s.mix=mix
        if mix: s.enc=nn.TransformerEncoder(nn.TransformerEncoderLayer(d,h,dim_feedforward=2*d,dropout=drop,batch_first=True,norm_first=True),2)
        else: s.enc=nn.Sequential(nn.Linear(d,2*d),nn.GELU(),nn.Dropout(drop),nn.Linear(2*d,d),nn.GELU(),nn.Linear(d,d))
        s.out=nn.Linear(d,1); nn.init.zeros_(s.out.weight); nn.init.zeros_(s.out.bias); s.bias=nn.Parameter(torch.zeros(1))
    def forward(s,x):
        B=x.shape[0]; xs=x.reshape(B*NS,G,M)
        e=s.fe(St)[None]+xs[...,None]*s.fd(St)[None]+s.le.weight[None,:,None,:]; e=e.reshape(B*NS*G,M,-1)
        if s.mix:
            km=Mt[None].expand(B*NS,-1,-1).reshape(B*NS*G,M); CH=1024; parts=[]
            for j in range(0,e.shape[0],CH):
                if s.training: parts.append(torch.utils.checkpoint.checkpoint(lambda a,b_: s.enc(a,src_key_padding_mask=b_),e[j:j+CH],km[j:j+CH],use_reentrant=False))
                else: parts.append(s.enc(e[j:j+CH],src_key_padding_mask=km[j:j+CH]))
            hh=torch.cat(parts,0)
        else: hh=s.enc(e)
        o=s.out(hh).reshape(B*NS,G,M).masked_fill(Mt[None],0.0)
        return (o*xs).sum((-1,-2)).reshape(B,NS).sum(-1)/SCALE+s.bias
def batch_x(ib):
    if ARM=="P5h":
        h=torch.tensor(H[ib][:,IDX,:].astype(np.float32))
        h=(h-Pt[IDX][None,:,None])/Sf[IDX][None,:,None]
        return h.permute(0,2,1).reshape(len(ib),2,G,M).to(dev)
    return torch.tensor(D[ib][:,IDX].astype(np.float32).reshape(len(ib),1,G,M),device=dev)
def predict(net,ids,bs=32):
    net.eval(); out=[]
    with torch.no_grad(), torch.autocast("cuda",dtype=torch.float16):
        for j in range(0,len(ids),bs): out.append(net(batch_x(ids[j:j+bs])).float().cpu().numpy())
    net.train(); return np.concatenate(out)
itr=np.where(tr)[0]; iva=np.where(va)[0]; ite=np.where(te)[0]; sub=np.random.default_rng(7).choice(itr,5000,replace=False)
LR,DROP,WD,BS,MAXEP,PATL=3e-5,0.1,1e-2,16,4,4
net=Net(drop=DROP,mix=(ARM!="P6")).to(dev); opt=torch.optim.AdamW(net.parameters(),lr=LR,weight_decay=WD); scaler=torch.cuda.amp.GradScaler()
p1_sub=r2(r[sub],base_pred[sub]); v0=r2(r[va],base_pred[va]+predict(net,iva)); print(f"{ARM} step0 valid_r2={v0:.6f} P1_train5k={p1_sub:.6f}",flush=True)
def step(ib):
    xb=batch_x(ib); yb=torch.tensor(rr[ib],device=dev)
    with torch.autocast("cuda",dtype=torch.float16): loss=((net(xb).float()-yb)**2).mean()
    opt.zero_grad(); scaler.scale(loss).backward(); scaler.unscale_(opt); torch.nn.utils.clip_grad_norm_(net.parameters(),1.0); scaler.step(opt); scaler.update()
net.train()
if MODE in ("smoke","bench"):
    perm=np.random.default_rng(1).permutation(itr); t0=time.time(); nb=200 if MODE=="smoke" else 30
    for b in range(nb): step(np.sort(perm[b*BS:(b+1)*BS]))
    el=(time.time()-t0)/nb; vt=r2(r[sub],base_pred[sub]+predict(net,sub)); v=r2(r[va],base_pred[va]+predict(net,iva))
    res=dict(arm=ARM,mode=MODE,s_per_batch=el,est_min_per_epoch=el*(len(itr)//BS)/60,train5k=vt,P1_train5k=p1_sub,valid=v,P1_valid=v0,gpu_gb=torch.cuda.max_memory_allocated()/1e9,
             smoke_pass=bool(vt>=p1_sub-0.005 and v>=v0-0.005))
    json.dump(res,open(O+f"/{MODE}_seed_{ARM}.json","w"),indent=1); print(json.dumps(res),flush=True); sys.exit(0 if (MODE=="bench" or res["smoke_pass"]) else 7)
bestv=v0; pat=0; hist=[]; stop=False; torch.save(net.state_dict(),P+f"/att_seed_{ARM}_best.pt")
for ep in range(MAXEP):
    perm=np.random.default_rng(20260927+ep).permutation(itr); nb=len(perm)//BS; CHECK=max(1,nb//4); t0=time.time()
    for b in range(nb):
        step(np.sort(perm[b*BS:(b+1)*BS]))
        if (b+1)%CHECK==0:
            v=r2(r[va],base_pred[va]+predict(net,iva)); vt=r2(r[sub],base_pred[sub]+predict(net,sub)); hist.append((ep,b+1,round(v,6),round(vt,6)))
            print(f"{ARM} ep={ep} step={b+1}/{nb} valid_r2={v:.6f} train5k_r2={vt:.6f} {time.time()-t0:.0f}s",flush=True)
            if v>bestv+1e-6: bestv=v; pat=0; torch.save(net.state_dict(),P+f"/att_seed_{ARM}_best.pt")
            else: pat+=1
            if pat>=PATL: stop=True; break
    if stop: break
net.load_state_dict(torch.load(P+f"/att_seed_{ARM}_best.pt"))
pred=base_pred.copy().astype(np.float32); pred[iva]+=predict(net,iva); pred[ite]+=predict(net,ite); np.save(P+f"/pred_att_seed_{ARM}.npy",pred)
json.dump(dict(arm=ARM,cap=CAP,genes=G,U=len(U),fs=fs,step0_valid=v0,best_valid_r2=bestv,valid_gain=bestv-v0,hist=hist,
               config=dict(lr=LR,drop=DROP,wd=WD,bs=BS,maxep=MAXEP,patience_checks=PATL,clip=1.0,scale=SCALE)),open(O+f"/train_att_seed_{ARM}.json","w"),indent=1)
print("DONE",flush=True)
