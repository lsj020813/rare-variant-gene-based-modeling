
import torch, torch.nn as nn
print(torch.__version__)
torch.manual_seed(0)
d=32; M=1024; B=6
enc=nn.TransformerEncoder(nn.TransformerEncoderLayer(d,4,dim_feedforward=2*d,dropout=0.1,batch_first=True,norm_first=True),2)
x=torch.randn(B,M,d)
km=torch.zeros(B,M,dtype=torch.bool); km[:,700:]=True
enc.train()
for m in enc.modules():
    if isinstance(m,nn.Dropout): m.p=0.0
for l in enc.layers: l.dropout.p=0; l.dropout1.p=0; l.dropout2.p=0; l.self_attn.dropout=0.0
ya=enc(x,src_key_padding_mask=km)
enc.eval()
with torch.no_grad(): yb=enc(x,src_key_padding_mask=km)
v=~km
print("max|train-eval| valid tokens", float((ya-yb)[v].abs().max()), "scale", float(ya[v].abs().mean()))
print("eval pad tokens mean|.|", float(yb[km].abs().mean()), "train pad mean", float(ya[km].abs().mean()))

with torch.no_grad(), torch.autocast("cpu",dtype=torch.bfloat16): yc=enc(x,src_key_padding_mask=km)
print("max|train-eval_autocast| valid", float((ya-yc.float())[v].abs().max()))
