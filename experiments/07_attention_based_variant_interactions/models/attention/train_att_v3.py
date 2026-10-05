#!/usr/bin/env python
import sys, json, time, itertools, numpy as np, torch, torch.nn as nn, torch.utils.checkpoint
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "data_prep"))
from prs_common import *
TAG=sys.argv[1]; ARM=sys.argv[2]; BENCH=len(sys.argv)>3
torch.manual_seed(20260926); np.random.seed(20260926)
r,prs,sp=base(TAG); tr=sp=="train"; va=sp=="valid"; te=sp=="test"
A=np.column_stack([np.ones(len(r)),prs]); ba=np.linalg.lstsq(A[tr],r[tr],rcond=None)[0]; base_pred=A@ba; rr=(r-base_pred).astype(np.float32)
Ls=loci(); F,cols=feats(); M=1024; nL=len(Ls); n=len(r)
X=np.zeros((n,nL,M),dtype=np.float16); S=None; mask=np.ones((nL,M),dtype=bool)
stat=[]
for i,L in enumerate(Ls):
    m=meta(L); D=dosage(L,tr); k=D.shape[1]; X[:,i,:k]=D.astype(np.float16); mask[i,:k]=False
    keys=[x["key"] for x in m]
    ann,names=encode_feats(F,cols,keys) if F else (np.zeros((k,0)),[])
    extra=np.column_stack([np.log10(np.array([float(x["maf_train"]) for x in m])),[float(x["r2"]) for x in m],[float(x["typed"]) for x in m],[float(x["relpos"]) for x in m]])
    stat.append((ann,names,extra))
allnames=sorted(set(nm for a,names,e in stat for nm in names))
fs=len(allnames)+4; S=np.zeros((nL,M,fs),dtype=np.float32)
for i,(ann,names,extra) in enumerate(stat):
    k=extra.shape[0]; idx=[allnames.index(nm) for nm in names]
    if idx:
        Si=S[i]; Si[:k, idx]=ann
    S[i,:k,-4:]=extra
valid_tok=~mask
mu=S[valid_tok].mean(0); sd=S[valid_tok].std(0); sd[sd<1e-6]=1; S=(S-mu)/sd; S[mask]=0
dev=torch.device("cuda")
St=torch.tensor(S,device=dev); Mt=torch.tensor(mask,device=dev)
class Net(nn.Module):
    def __init__(s, d=32, h=4, drop=0.0, mix=True):
        super().__init__()
        s.fe=nn.Sequential(nn.Linear(fs,d),nn.GELU(),nn.Linear(d,d)); s.fd=nn.Linear(fs,d); s.le=nn.Embedding(nL,d)
        s.mix=mix
        if mix:
            s.enc=nn.TransformerEncoder(nn.TransformerEncoderLayer(d,h,dim_feedforward=2*d,dropout=drop,batch_first=True,norm_first=True),2)
        else:
            s.enc=nn.Sequential(nn.Linear(d,2*d),nn.GELU(),nn.Dropout(drop),nn.Linear(2*d,d),nn.GELU(),nn.Linear(d,d))
        s.out=nn.Linear(d,1); nn.init.zeros_(s.out.weight); nn.init.zeros_(s.out.bias)
    def forward(s, x):
        B=x.shape[0]
        e=s.fe(St)[None]+x[...,None]*s.fd(St)[None]+s.le.weight[None,:,None,:]
        e=e.reshape(B*nL,M,-1)
        if s.mix:
            km=Mt[None].expand(B,-1,-1).reshape(B*nL,M); CH=42; parts=[]
            for j in range(0,B*nL,CH):
                if s.training:
                    parts.append(torch.utils.checkpoint.checkpoint(lambda a,b: s.enc(a,src_key_padding_mask=b), e[j:j+CH], km[j:j+CH], use_reentrant=False))
                else:
                    parts.append(s.enc(e[j:j+CH], src_key_padding_mask=km[j:j+CH]))
            hh=torch.cat(parts,0)
        else: hh=s.enc(e)
        o=s.out(hh).reshape(B,nL,M)
        o=o.masked_fill(Mt[None],0.0)
        return (o*x).sum(-1).sum(-1)/np.sqrt(M)
def predict(net, idx, bs=64):
    net.eval(); out=[]
    with torch.no_grad(), torch.autocast("cuda",dtype=torch.float16):
        for j in range(0,len(idx),bs):
            xb=torch.tensor(X[idx[j:j+bs]].astype(np.float32),device=dev); out.append(net(xb).float().cpu().numpy())
    return np.concatenate(out)
itr=np.where(tr)[0]; iva=np.where(va)[0]; ite=np.where(te)[0]
BS=16 if ARM=="P5" else 64
grid=list(itertools.product([3e-4,1e-3],[0.0,0.1],[0.0,1e-4]))
if BENCH: grid=grid[:1]
log={"tag":TAG,"arm":ARM,"n_loci":nL,"fs":fs,"grid":[]}
best=None
for lr,drop,wd in grid:
    torch.manual_seed(20260926)
    net=Net(drop=drop,mix=(ARM=="P5")).to(dev); opt=torch.optim.AdamW(net.parameters(),lr=lr,weight_decay=wd); scaler=torch.cuda.amp.GradScaler()
    bestv=-1e9; pat=0; hist=[]
    for ep in range(30):
        net.train(); perm=np.random.default_rng(20260926+ep).permutation(itr); t0=time.time()
        nb=len(perm)//BS
        for b in range(nb):
            ib=np.sort(perm[b*BS:(b+1)*BS])
            xb=torch.tensor(X[ib].astype(np.float32),device=dev); yb=torch.tensor(rr[ib],device=dev)
            with torch.autocast("cuda",dtype=torch.float16): loss=((net(xb).float()-yb)**2).mean()
            opt.zero_grad(); scaler.scale(loss).backward(); scaler.step(opt); scaler.update()
            if BENCH and b==49:
                el=time.time()-t0; json.dump({"bench_s_per_batch":el/50,"batches_per_epoch":nb,"est_min_per_epoch":el/50*nb/60,"gpu_mem_gb":torch.cuda.max_memory_allocated()/1e9,"fs":fs,"n_loci":nL},open(f"{W}/out/bench_{ARM}.json","w")); print("bench",el/50*nb/60); sys.exit(0)
        pv=predict(net,iva); v=r2(r[va],base_pred[va]+pv); hist.append(round(v,6))
        print(f"{ARM} lr={lr} drop={drop} wd={wd} ep={ep} valid_r2={v:.6f} {time.time()-t0:.0f}s",flush=True)
        if v>bestv+1e-6: bestv=v; pat=0; torch.save(net.state_dict(),f"{P}/model/{TAG}_{ARM}_cur.pt")
        else:
            pat+=1
            if pat>=3: break
    log["grid"].append({"lr":lr,"drop":drop,"wd":wd,"best_valid_r2":bestv,"epochs":len(hist),"hist":hist})
    if best is None or bestv>best[0]:
        best=(bestv,(lr,drop,wd)); import shutil; shutil.copy(f"{P}/model/{TAG}_{ARM}_cur.pt",f"{P}/model/{TAG}_{ARM}_best.pt")
net=Net(drop=best[1][1],mix=(ARM=="P5")).to(dev); net.load_state_dict(torch.load(f"{P}/model/{TAG}_{ARM}_best.pt"))
pred=np.zeros(n,dtype=np.float32); pred[iva]=predict(net,iva); pred[ite]=predict(net,ite)
np.save(f"{P}/model/pred_{TAG}_{ARM}.npy", (base_pred+pred).astype(np.float32))
log["selected"]={"lr":best[1][0],"drop":best[1][1],"wd":best[1][2],"valid_r2":best[0]}
json.dump(log,open(f"{W}/out/train_{TAG}_{ARM}.json","w"),indent=1); print(json.dumps(log["selected"]))
