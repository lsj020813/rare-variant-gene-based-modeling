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
N_GENES = _required_number("N_GENES", int, True)
import sys, time, resource, torch, torch.nn as nn, torch.utils.checkpoint, json
dev=sys.argv[1]; th=int(sys.argv[2]); torch.set_num_threads(th); torch.manual_seed(0)
G,M,fs,B=N_GENES,128,101,16
St=torch.randn(G,M,fs,device=dev); Mt=torch.zeros(G,M,dtype=torch.bool,device=dev); Mt[:,110:]=True
class Net(nn.Module):
    def __init__(s,NS,d=16,h=2,drop=0.2):
        super().__init__(); s.NS=NS
        s.fe=nn.Sequential(nn.Linear(fs,d),nn.GELU(),nn.Linear(d,d)); s.fd=nn.Linear(fs,d); s.le=nn.Embedding(G,d)
        s.enc=nn.TransformerEncoder(nn.TransformerEncoderLayer(d,h,dim_feedforward=2*d,dropout=drop,batch_first=True,norm_first=True),1)
        s.out=nn.Linear(d,1)
    def forward(s,x):
        Bn=x.shape[0]; xs=x.reshape(Bn*s.NS,G,M)
        e=s.fe(St)[None]+xs[...,None]*s.fd(St)[None]+s.le.weight[None,:,None,:]; e=e.reshape(Bn*s.NS*G,M,-1)
        km=Mt[None].expand(Bn*s.NS,-1,-1).reshape(Bn*s.NS*G,M); parts=[]
        for j in range(0,e.shape[0],1024):
            parts.append(torch.utils.checkpoint.checkpoint(lambda a,b_: s.enc(a,src_key_padding_mask=b_),e[j:j+1024],km[j:j+1024],use_reentrant=False))
        o=s.out(torch.cat(parts,0)).reshape(Bn*s.NS,G,M).masked_fill(Mt[None],0.0)
        return (o*xs).sum((-1,-2)).reshape(Bn,s.NS).sum(-1)
res={}
for NS in (1,2):
    net=Net(NS).to(dev); opt=torch.optim.AdamW(net.parameters(),lr=3e-5); x=torch.randn(B,NS,G,M,device=dev); y=torch.randn(B,device=dev)
    ts=[]
    for it in range(10):
        if dev=="cuda": torch.cuda.synchronize()
        t=time.time()
        if dev=="cuda":
            with torch.autocast("cuda",dtype=torch.float16): loss=((net(x).float()-y)**2).mean()
        else: loss=((net(x)-y)**2).mean()
        opt.zero_grad(); loss.backward(); opt.step()
        if dev=="cuda": torch.cuda.synchronize()
        ts.append(time.time()-t)
    ts=sorted(ts[2:]); res[f"NS{NS}_s_per_batch"]=round(ts[len(ts)//2],3)
print(json.dumps(dict(dev=dev,threads=th,peak_rss_gb=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6,2),**res)))
