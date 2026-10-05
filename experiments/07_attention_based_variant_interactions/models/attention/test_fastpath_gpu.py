
import torch, torch.nn as nn
torch.manual_seed(0); dev="cuda"
d=32; M=1024; B=6
enc=nn.TransformerEncoder(nn.TransformerEncoderLayer(d,4,dim_feedforward=2*d,dropout=0.0,batch_first=True,norm_first=True),2).to(dev)
x=torch.randn(B,M,d,device=dev); km=torch.zeros(B,M,dtype=torch.bool,device=dev); km[:,700:]=True; v=~km
enc.train()
with torch.autocast("cuda",dtype=torch.float16): ya=enc(x,src_key_padding_mask=km).float()
enc.eval()
with torch.no_grad(): yb=enc(x,src_key_padding_mask=km).float()
with torch.no_grad(), torch.autocast("cuda",dtype=torch.float16): yc=enc(x,src_key_padding_mask=km).float()
print("train16 vs eval32", float((ya-yb)[v].abs().max()), "| train16 vs eval16", float((ya-yc)[v].abs().max()), "| scale", float(ya[v].abs().mean()))
print("nan eval16", int(torch.isnan(yc).sum()), "pad eval16 mean", float(yc[km].abs().mean()))
