#!/usr/bin/env python
import os as _os
import math as _number_math
def _required_number(name, cast, positive=False):
    raw = _os.environ.get(name, "")
    if not raw.strip():
        raise ValueError(name + " must be set and nonblank")
    try:
        value = cast(raw)
    except (ValueError, OverflowError):
        raise ValueError(name + " has an invalid numeric value") from None
    if isinstance(value, float) and not _number_math.isfinite(value):
        raise ValueError(name + " must be finite")
    if positive and value <= 0:
        raise ValueError(name + " must be positive")
    return value
if not _os.environ.get("PROJECT_ROOT", "").strip():
    raise ValueError("PROJECT_ROOT must be set and nonblank")
EXPECTED_BASELINE_R2 = _required_number("EXPECTED_BASELINE_R2", float, False)
import sys, os, json, time, datetime, csv, gzip, glob, collections, numpy as np, torch, torch.nn as nn, torch.utils.checkpoint
W=(_os.environ["PROJECT_ROOT"] + '/work/prs'); P0=W+"/private/hei"; P=W+"/private/hei_seed"; O=W+"/out/hei_seed"; SEL=W+"/out/hei/sel_seed"
ARM=sys.argv[1]; assert ARM in ("P6res","P5res","P5hres"); BA=ARM[:-3]; LOGTEST=os.environ.get("HEI_LOGTEST")=="1"; MODE=sys.argv[2] if len(sys.argv)>2 else "train"
torch.manual_seed(20260927); np.random.seed(20260927)
d=np.load(P0+"/base_hei.npz",allow_pickle=True); r=d["r"].astype(np.float64); prs=d["prs"].astype(np.float64); sp=d["split"]
tr=sp=="train"; va=sp=="valid"; te=sp=="test"; n=len(r)
base_pred=np.load(P+"/p3ann128_oof.npy").astype(np.float64); assert base_pred.shape==(n,); rr=(r-base_pred).astype(np.float32)
def r2(y,p): return float(1-((y-p)**2).sum()/((y-y.mean())**2).sum())
from concurrent.futures import ThreadPoolExecutor
import resource
def now(): return datetime.datetime.now().strftime("%m-%d %H:%M:%S")
def rss(): return round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6,1)
T_LOAD=time.time(); print(f"[{now()}] {ARM} loading start",flush=True)
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
if BA=="P5h":
    loc={}
    for f in glob.glob(O+"/gt/k_chr*.tsv"):
        N=os.path.basename(f)[5:-4]
        for j,row in enumerate(csv.DictReader(open(f),delimiter="\t")): loc[row["key"]]=(N,j)
    assert all(k in loc for k in U), f"gt missing {sum(k not in loc for k in U)}"
    H=np.zeros((n,len(U),2),np.int8); byc=collections.defaultdict(list)
    for k in U: byc[loc[k][0]].append((loc[k][1],ui[k]))
    def _lh(item):
        N,lst=item; t=time.time(); Hc=np.load(f"{P}/gt/h_chr{N}.npy"); cols=np.array([a for a,_ in lst]); dst=np.array([b for _,b in lst])
        H[:,dst,:]=Hc[:,cols,:]; del Hc; return N,len(cols),time.time()-t
    with ThreadPoolExecutor(8) as ex:
        for N,m_,dt in ex.map(_lh,sorted(byc.items(),key=lambda x:-len(x[1]))): print(f"[{now()}] load gt chr{N} cols={m_} {dt:.0f}s rss={rss()}GB",flush=True)
    if "22" in byc:
        lst=byc["22"]; rs_=np.random.default_rng(3).choice(len(lst),min(200,len(lst)),replace=False); Hm=np.load(f"{P}/gt/h_chr22.npy",mmap_mode="r")
        for q in rs_: a,b_=lst[q]; assert np.array_equal(np.asarray(Hm[:,a,:]),H[:,b_,:]), f"gt load mismatch chr22 col {a}"
        del Hm; print(f"[{now()}] gt load equality check chr22 OK ({len(rs_)} cols)",flush=True)
    pf=H[tr].mean(axis=(0,2)).astype(np.float32); sf=np.sqrt(np.maximum(pf*(1-pf),1e-4)).astype(np.float32)
    Pt=torch.tensor(pf); Sf=torch.tensor(sf)
else:
    loc={}
    for f in glob.glob(O+"/meta/meta_chr*.tsv"):
        N=os.path.basename(f)[8:-4]
        for j,row in enumerate(csv.DictReader(open(f),delimiter="\t")): loc[row["key"]]=(N,j)
    D=np.zeros((n,len(U)),np.float16); byc=collections.defaultdict(list)
    for k in U: byc[loc[k][0]].append((loc[k][1],ui[k]))
    def _std(X):
        X=X.astype(np.float32); mu=X[tr].mean(0); sd=X[tr].std(0); sd[sd<1e-6]=1; return ((X-mu)/sd).astype(np.float16)
    def _ld(item):
        N,lst=item; t=time.time(); Gm=np.load(f"{P}/geno_chr{N}.npy"); cols=np.array([a for a,_ in lst]); dst=np.array([b for _,b in lst])
        for s in range(0,len(cols),4096): D[:,dst[s:s+4096]]=_std(Gm[:,cols[s:s+4096]])
        del Gm; return N,len(cols),time.time()-t
    with ThreadPoolExecutor(8) as ex:
        for N,m_,dt in ex.map(_ld,sorted(byc.items(),key=lambda x:-len(x[1]))): print(f"[{now()}] load dosage chr{N} cols={m_} {dt:.0f}s rss={rss()}GB",flush=True)
    if "22" in byc:
        lst=byc["22"]; rs_=np.random.default_rng(3).choice(len(lst),min(200,len(lst)),replace=False); Gm=np.load(f"{P}/geno_chr22.npy",mmap_mode="r")
        for q in rs_: a,b_=lst[q]; assert np.array_equal(_std(np.asarray(Gm[:,[a]]))[:,0],D[:,b_]), f"dosage load mismatch chr22 col {a}"
        del Gm; print(f"[{now()}] dosage load equality check chr22 OK ({len(rs_)} cols)",flush=True)
print(f"[{now()}] genotype stage done {time.time()-T_LOAD:.0f}s rss={rss()}GB",flush=True)
F={}; cols_f=None
for f in sorted(glob.glob(O+"/feat/feat_chr*.tsv.gz")):
    with gzip.open(f,"rt") as g:
        h=g.readline().rstrip("\n").split("\t")
        if cols_f is None: cols_f=[c for c in h if c not in ("key","chr","pos38","domain") and "re2g" not in c.lower() and not c.endswith("_target")]
        ix=[h.index(c) for c in cols_f]
        for l in g:
            p_=l.rstrip("\n").split("\t"); F[p_[0]]=[p_[i] for i in ix]
print(f"[{now()}] feature files read {time.time()-T_LOAD:.0f}s",flush=True); fmiss=sum(k not in F for k in U); assert fmiss==0, f"features missing for {fmiss} tokens"
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
NS=2 if BA=="P5h" else 1; SCALE=float(np.sqrt(G*M*NS))
print(f"[{now()}] data loaded total {time.time()-T_LOAD:.0f}s rss={rss()}GB",flush=True); print(json.dumps(dict(arm=ARM,mode=MODE,genes=G,M=M,U=len(U),fs=fs,n=n,strands=NS,scale=SCALE)),flush=True)
class Net(nn.Module):
    def __init__(s,d=16,h=2,drop=0.2,mix=True):
        super().__init__()
        s.fe=nn.Sequential(nn.Linear(fs,d),nn.GELU(),nn.Linear(d,d)); s.fd=nn.Linear(fs,d); s.le=nn.Embedding(G,d); s.mix=mix
        if mix: s.enc=nn.TransformerEncoder(nn.TransformerEncoderLayer(d,h,dim_feedforward=2*d,dropout=drop,batch_first=True,norm_first=True),1)
        else: s.enc=nn.Sequential(nn.Linear(d,2*d),nn.GELU(),nn.Dropout(drop),nn.Linear(2*d,d))
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
    if BA=="P5h":
        h=torch.tensor(H[ib][:,IDX,:].astype(np.float32))
        h=(h-Pt[IDX][None,:,None])/Sf[IDX][None,:,None]
        return h.permute(0,2,1).reshape(len(ib),2,G,M).to(dev)
    return torch.tensor(D[ib][:,IDX].astype(np.float32).reshape(len(ib),1,G,M),device=dev)
def gpu(): return round(torch.cuda.memory_allocated()/1e9,2) if torch.cuda.is_available() else 0
STATUS=O+f"/status_res_{ARM}"+("_logtest" if LOGTEST else "")+".json"
def status(**kw):
    kw.update(arm=ARM,time=now()); json.dump(kw,open(STATUS+".tmp","w"),indent=1); os.replace(STATUS+".tmp",STATUS)
def predict(net,ids,bs=32,tag=""):
    net.eval(); out=[]; t=time.time(); nbp=(len(ids)+bs-1)//bs
    with torch.no_grad(), torch.autocast("cuda",dtype=torch.float16):
        for k,j in enumerate(range(0,len(ids),bs)):
            out.append(net(batch_x(ids[j:j+bs])).float().cpu().numpy())
            if tag and (k+1)%(5 if LOGTEST else 100)==0:
                el=time.time()-t; print(f"[{now()}] {ARM} predict {tag} {k+1}/{nbp} ({100*(k+1)/nbp:.0f}%) eta {el/(k+1)*(nbp-k-1)/60:.1f} min",flush=True)
    if tag: print(f"[{now()}] {ARM} predict {tag} done n={len(ids)} {time.time()-t:.0f}s",flush=True)
    net.train(); return np.concatenate(out)
itr=np.where(tr)[0]; iva=np.where(va)[0]; ite=np.where(te)[0]; sub=np.random.default_rng(7).choice(itr,10000,replace=False)
if LOGTEST: iva=iva[:640]; ite=ite[:320]; sub=sub[:320]; MAXEP_OVERRIDE=1
LR,DROP,WD,BS,MAXEP,PATL=3e-5,0.2,1e-2,16,4,4
if LOGTEST: MAXEP=1
net=Net(drop=DROP,mix=(BA!="P6")).to(dev); opt=torch.optim.AdamW(net.parameters(),lr=LR,weight_decay=WD); scaler=torch.cuda.amp.GradScaler()
p1_sub=r2(r[sub],base_pred[sub]); v0=r2(r[va],base_pred[va]+predict(net,iva)); print(f"{ARM} step0 valid_r2={v0:.6f} P3oof_train10k={p1_sub:.6f}",flush=True); assert LOGTEST or abs(v0-EXPECTED_BASELINE_R2)<3e-4, v0
status(stage="start",step0_valid=v0,step0_train10k=p1_sub,mode=MODE)
def step(ib):
    xb=batch_x(ib); yb=torch.tensor(rr[ib],device=dev)
    with torch.autocast("cuda",dtype=torch.float16): loss=((net(xb).float()-yb)**2).mean()
    opt.zero_grad(); scaler.scale(loss).backward(); scaler.unscale_(opt); torch.nn.utils.clip_grad_norm_(net.parameters(),1.0); scaler.step(opt); scaler.update()
    return float(loss.detach())
net.train()
if MODE in ("smoke","bench"):
    perm=np.random.default_rng(1).permutation(itr); t0=time.time(); nb=200 if MODE=="smoke" else 30
    for b in range(nb):
        l=step(np.sort(perm[b*BS:(b+1)*BS]))
        if (b+1)%25==0: print(f"[{now()}] {ARM} smoke batch {b+1}/{nb} loss={l:.6f} {(time.time()-t0)/(b+1):.2f}s/batch gpu={gpu()}GB",flush=True)
    el=(time.time()-t0)/nb; vt=r2(r[sub],base_pred[sub]+predict(net,sub)); v=r2(r[va],base_pred[va]+predict(net,iva))
    res=dict(arm=ARM,mode=MODE,s_per_batch=el,est_min_per_epoch=el*(len(itr)//BS)/60,train5k=vt,P1_train5k=p1_sub,valid=v,P1_valid=v0,gpu_gb=torch.cuda.max_memory_allocated()/1e9,
             smoke_pass=bool(vt>=p1_sub-0.005 and v>=v0-0.005)); status(stage="smoke_done",**{k:(float(x) if isinstance(x,(float,np.floating)) else x) for k,x in res.items()})
    json.dump(res,open(O+f"/{MODE}_res_{ARM}.json","w"),indent=1); print(json.dumps(res),flush=True); sys.exit(0 if (MODE=="bench" or res["smoke_pass"]) else 7)
T_START=time.time(); LH=[]; bestv=v0; pat=0; hist=[]; stop=False; torch.save(net.state_dict(),P+f"/att_res_{ARM}"+("_logtest" if LOGTEST else "")+"_best.pt")
for ep in range(MAXEP):
    perm=np.random.default_rng(20260927+ep).permutation(itr); nb=len(perm)//BS; CHECK=max(1,nb//4); t0=time.time()
    if LOGTEST: nb=130; CHECK=60; tb=time.time(); bstart=0
    print(f"[{now()}] {ARM} === epoch {ep} start, {nb} batches, check every {CHECK} ===",flush=True)
    for b in range(nb):
        lsum=locals().get("lsum",0.0)+step(np.sort(perm[b*BS:(b+1)*BS]))
        if (b+1)%50==0:
            el=time.time()-tb; spb=el/((b+1)-bstart); nxt=CHECK-((b+1)%CHECK) if (b+1)%CHECK else 0
            print(f"[{now()}] {ARM} ep={ep} batch={b+1}/{nb} ({100*(b+1)/nb:.1f}% of epoch) loss50={lsum/50:.6f} {spb:.2f}s/batch next_check_in={nxt} batches (~{nxt*spb/60:.0f} min + eval) gpu={gpu()}GB best={bestv:.6f} patience={pat}/{PATL}",flush=True)
            LH.append(lsum/50)
            if time.time()-globals().get("last_sum",T_START)>=(20 if LOGTEST else 1800):
                last_sum=time.time(); k=len(LH)//2
                trend=(sum(LH[k:])/max(1,len(LH[k:])))-(sum(LH[:k])/max(1,k)) if len(LH)>=4 else float("nan")
                lastc=(f"last check valid={hist[-1][2]:.6f} (P3 {v0:.6f}, {hist[-1][2]-v0:+.6f})" if hist else "no check yet")
                rem_ep=(nb-(b+1))*spb/3600
                print(f"[{now()}] {ARM} ##### 30-MIN SUMMARY: ep {ep}/{MAXEP} batch {b+1}/{nb} ({100*(b+1)/nb:.1f}%), elapsed {(time.time()-T_START)/3600:.2f} h, loss mean(last half)-mean(first half) since last summary = {trend:+.6f}, {lastc}, best={bestv:.6f}, patience={pat}/{PATL}, next check ~{nxt*spb/60:.0f} min + eval, epoch end ~{rem_ep:.1f} h + evals #####",flush=True); LH=[]
            status(stage="train",epoch=ep,batch=b+1,nb=nb,pct_epoch=round(100*(b+1)/nb,1),loss50=lsum/50,s_per_batch=spb,next_check_in_batches=nxt,best_valid=bestv,patience=pat,hist=hist); lsum=0.0
        if (b+1)%CHECK==0:
            print(f"[{now()}] {ARM} ep={ep} batch={b+1} === CHECK {len(hist)+1}: validating (valid {len(iva)} + train10k) ===",flush=True); status(stage="validating",epoch=ep,batch=b+1,nb=nb,best_valid=bestv,patience=pat,hist=hist)
            tv=time.time(); v=r2(r[va],base_pred[va]+predict(net,iva,tag="valid")); vt=r2(r[sub],base_pred[sub]+predict(net,sub,tag="train10k")); hist.append((ep,b+1,round(v,6),round(vt,6)))
            print(f"{ARM} ep={ep} step={b+1}/{nb} valid_r2={v:.6f} train5k_r2={vt:.6f} {time.time()-t0:.0f}s",flush=True)
            improved=v>bestv+1e-6
            if improved: bestv=v; pat=0; torch.save(net.state_dict(),P+f"/att_res_{ARM}"+("_logtest" if LOGTEST else "")+"_best.pt")
            else: pat+=1
            torch.save(net.state_dict(),P+f"/att_res_{ARM}"+("_logtest" if LOGTEST else "")+f"_ep{ep}_b{b+1}.pt")
            print(f"[{now()}] {ARM} >>> CHECK {len(hist)} RESULT ep={ep} batch={b+1}: valid={v:.6f} (vs P3 {v0:.6f}: {v-v0:+.6f}) train10k={vt:.6f} (vs P3oof {p1_sub:.6f}: {vt-p1_sub:+.6f}) best={bestv:.6f} {'IMPROVED' if improved else 'no improvement'} patience={pat}/{PATL} eval_took={time.time()-tv:.0f}s",flush=True)
            status(stage="train",epoch=ep,batch=b+1,nb=nb,last_valid=v,last_train10k=vt,step0_valid=v0,step0_train10k=p1_sub,best_valid=bestv,patience=pat,hist=hist)
            if pat>=PATL: stop=True; break
    if stop: break
net.load_state_dict(torch.load(P+f"/att_res_{ARM}"+("_logtest" if LOGTEST else "")+"_best.pt"))
print(f"[{now()}] {ARM} === training stopped (stop={stop}); final predictions with best model ===",flush=True); status(stage="final_predict",best_valid=bestv,hist=hist)
pred=base_pred.copy().astype(np.float32); pred[iva]+=predict(net,iva,tag="final_valid"); pred[ite]+=predict(net,ite,tag="final_test(saved only, not scored)"); np.save(P+f"/pred_att_res_{ARM}"+("_logtest" if LOGTEST else "")+".npy",pred)
json.dump(dict(arm=ARM,cap=CAP,genes=G,U=len(U),fs=fs,step0_valid=v0,best_valid_r2=bestv,valid_gain=bestv-v0,hist=hist,
               config=dict(lr=LR,drop=DROP,wd=WD,bs=BS,maxep=MAXEP,patience_checks=PATL,clip=1.0,scale=SCALE)),open(O+f"/train_att_res_{ARM}"+("_logtest" if LOGTEST else "")+".json","w"),indent=1)
status(stage="done",best_valid=bestv,valid_gain=bestv-v0,hist=hist); print(f"[{now()}] DONE",flush=True)
