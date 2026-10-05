#!/usr/bin/env python3
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import os, sys, time, json, math, argparse, hashlib
os.environ.setdefault("HF_HOME",_config_path("${PROJECT_ROOT}/work/run_trackB/hf_cache"))
os.environ.setdefault("TMPDIR",_config_path("${PROJECT_ROOT}/work/tmp"))
import numpy as np, torch
from pyfaidx import Fasta
from enformer_pytorch import from_pretrained, SEQUENCE_LENGTH, TARGET_LENGTH
ap=argparse.ArgumentParser()
ap.add_argument("--variants", default=_config_path("${PROJECT_ROOT}/work/run_trackB/variants.tsv"))
ap.add_argument("--out", required=True)
ap.add_argument("--limit", type=int, default=0)
ap.add_argument("--fasta", default=_config_path("${PROJECT_ROOT}/work/run_trackB/ref/chr19.hg19.fa"))
ap.add_argument("--targets", default=_config_path("${PROJECT_ROOT}/work/run_trackB/ref/targets_human.txt"))
ap.add_argument("--manifest", default=_config_path("${PROJECT_ROOT}/work/run_trackB/score_manifest.json"))
a=ap.parse_args()
torch.set_num_threads(4)
L=SEQUENCE_LENGTH; NB=TARGET_LENGTH; BIN=128
assert L==196608 and NB==896
OFF=(L-NB*BIN)//2
CB=(L//2-OFF)//BIN
assert CB==448
C3=[CB-1,CB,CB+1]
KW=["liver","hepg2","hepatocyte","blood","k562","gm12878"]
desc=[]
with open(a.targets) as fh:
    hdr=fh.readline().rstrip("\n").split("\t"); di=hdr.index("description")
    for line in fh: desc.append(line.rstrip("\n").split("\t")[di])
assert len(desc)==5313
is_cage=np.array([d.upper().startswith("CAGE:") for d in desc]); is_dnase=np.array([d.upper().startswith("DNASE:") for d in desc])
kw=np.array([any(k in d.lower() for k in KW) for d in desc])
tis=(is_cage|is_dnase)&kw; tis_cage=is_cage&kw
man={"model":"EleutherAI/enformer-official-rough (enformer-pytorch 0.8.12, use_tf_gamma default for official weights)",
     "input_bp":L,"output_bins":NB,"bin_bp":BIN,"centre_bin":CB,"centre_bins_pm1":C3,"fasta":a.fasta,
     "fasta_sha256":hashlib.sha256(open(a.fasta,"rb").read()).hexdigest(),"targets_sha256":hashlib.sha256(open(a.targets,"rb").read()).hexdigest(),
     "keywords":KW,"n_tracks":5313,"n_cage":int(is_cage.sum()),"n_dnase":int(is_dnase.sum()),"n_tissue_subset":int(tis.sum()),"n_tissue_cage":int(tis_cage.sum()),
     "tissue_subset_idx":np.where(tis)[0].tolist(),
     "scores":{"S1_all_c3_log":"sum over all 5313 tracks, bins 447-449, |log1p(alt)-log1p(ref)|",
               "S2_all_full_log":"sum over all tracks, all 896 bins, |dlog1p|",
               "S3_tis_c3_log":"sum over tissue-subset CAGE+DNase tracks, bins 447-449, |dlog1p|",
               "S4_tis_full_log":"sum over tissue-subset tracks, all bins, |dlog1p|",
               "S1r_all_c3_raw":"same as S1 on raw scale |alt-ref|","S2r_all_full_raw":"same as S2 on raw scale",
               "sign_cage_c3":"fraction of ALL CAGE tracks with (sum bins 447-449 dlog1p)>0; majority sign = sign(frac-0.5)",
               "cage_tis_c3_signed":"signed sum over tissue-subset CAGE tracks, bins 447-449, dlog1p",
               "ref_c3_log":"sum over all tracks bins 447-449 log1p(ref) (covariate)"},
     "indel_handling":"alt = ref window with REF replaced by ALT, then re-cut to L keeping variant start at index L/2; bins downstream shift by len(alt)-len(ref) bp",
     "precision":"fp32, torch.no_grad, batch = [ref,alt] (2) with fallback to 1 on OOM",
     "torch":torch.__version__,"gpu":torch.cuda.get_device_name(0) if torch.cuda.is_available() else None}
if not os.path.exists(a.manifest):
    json.dump(man, open(a.manifest,"w"), indent=1)
print("tracks: all", 5313, "cage", man["n_cage"], "dnase", man["n_dnase"], "tissue_subset", man["n_tissue_subset"], "tissue_cage", man["n_tissue_cage"])
dev="cuda"
t0=time.time(); model=from_pretrained("EleutherAI/enformer-official-rough").to(dev).eval(); tload=time.time()-t0
print(f"model loaded {tload:.1f}s; params {sum(p.numel() for p in model.parameters())/1e6:.1f}M; gpu_alloc {torch.cuda.memory_allocated()/2**30:.2f}G")
fa=Fasta(a.fasta); chrom=[k for k in fa.keys()][0]; clen=len(fa[chrom]); print("fasta contig", chrom, "len", clen)
MAP={"A":0,"C":1,"G":2,"T":3}
def onehot(s):
    x=np.zeros((len(s),4),dtype=np.float32)
    for i,ch in enumerate(s):
        j=MAP.get(ch)
        if j is not None: x[i,j]=1.0
    return x
done=set()
if os.path.exists(a.out):
    with open(a.out) as fh:
        fh.readline()
        for line in fh: done.add(line.split("\t")[0])
    fo=open(a.out,"a")
else:
    fo=open(a.out,"w"); fo.write("key\tref_match\tn_frac\tt_sec\tS1_all_c3_log\tS2_all_full_log\tS3_tis_c3_log\tS4_tis_full_log\tS1r_all_c3_raw\tS2r_all_full_raw\tsign_cage_c3\tcage_tis_c3_signed\tref_c3_log\n"); fo.flush()
rows=[]
with open(a.variants) as fh:
    hdr=fh.readline().rstrip("\n").split("\t")
    for line in fh:
        t=dict(zip(hdr,line.rstrip("\n").split("\t"))); rows.append(t)
if a.limit: rows=rows[:a.limit]
todo=[r for r in rows if r["key"] not in done]
print("variants total", len(rows), "already done", len(done), "todo", len(todo))
PAD=64
def windows(pos,ref,alt):
    start0=pos-1-L//2
    s=start0-PAD; e=start0+L+PAD+len(alt)
    if s<0 or e>clen: return None
    reg=str(fa[chrom][s:e]).upper()
    i=PAD+L//2
    if reg[i:i+len(ref)]!=ref.upper(): return "MISMATCH"
    refw=reg[PAD:PAD+L]
    altreg=reg[:i]+alt.upper()+reg[i+len(ref):]
    altw=altreg[PAD:PAD+L]
    assert len(refw)==L and len(altw)==L
    return refw,altw
n_mm=0; n_oob=0; times=[]; peak=0
cage_idx=np.where(is_cage)[0]; tis_idx=np.where(tis)[0]; tisc_idx=np.where(tis_cage)[0]
def fwd(x):
    with torch.no_grad():
        return model(x)["human"].float().cpu().numpy()
for n,r in enumerate(todo):
    pos=int(r["pos"]); ref=r["ref"]; alt=r["alt"]; key=r["key"]
    w=windows(pos,ref,alt)
    if w is None: n_oob+=1; fo.write(f"{key}\tOOB\t\t\t\t\t\t\t\t\t\t\t\n"); fo.flush(); continue
    if w=="MISMATCH": n_mm+=1; fo.write(f"{key}\tMISMATCH\t\t\t\t\t\t\t\t\t\t\t\n"); fo.flush(); continue
    refw,altw=w; nfrac=refw.count("N")/L
    t1=time.time()
    x=torch.from_numpy(np.stack([onehot(refw),onehot(altw)])).to(dev)
    try:
        pred=fwd(x); pr,pa=pred[0],pred[1]
    except RuntimeError as ex:
        if "out of memory" not in str(ex): raise
        torch.cuda.empty_cache(); print("OOM at batch2 -> batch1"); pr=fwd(x[:1])[0]; pa=fwd(x[1:])[0]
    dt=time.time()-t1; times.append(dt); peak=max(peak, torch.cuda.max_memory_allocated())
    lr=np.log1p(pr); la=np.log1p(pa); dl=la-lr; dr=pa-pr
    dl3=dl[C3].sum(0)
    S1=float(np.abs(dl[C3]).sum()); S2=float(np.abs(dl).sum())
    S3=float(np.abs(dl[C3][:,tis_idx]).sum()); S4=float(np.abs(dl[:,tis_idx]).sum())
    S1r=float(np.abs(dr[C3]).sum()); S2r=float(np.abs(dr).sum())
    sgn=float((dl3[cage_idx]>0).mean()); ctis=float(dl3[tisc_idx].sum()); refc=float(lr[C3].sum())
    fo.write(f"{key}\tOK\t{nfrac:.4f}\t{dt:.3f}\t{S1:.6g}\t{S2:.6g}\t{S3:.6g}\t{S4:.6g}\t{S1r:.6g}\t{S2r:.6g}\t{sgn:.4f}\t{ctis:.6g}\t{refc:.6g}\n"); fo.flush()
    if (n+1)%25==0 or n<3:
        print(f"[{n+1}/{len(todo)}] t/var {np.mean(times):.2f}s (last {dt:.2f}) peak_gpu {peak/2**30:.2f}G mm {n_mm} oob {n_oob}", flush=True)
fo.close()
summ={"n_scored":len(times),"n_mismatch":n_mm,"n_oob":n_oob,"t_per_variant_mean_s":float(np.mean(times)) if times else None,
      "t_per_variant_median_s":float(np.median(times)) if times else None,"t_model_load_s":tload,"peak_gpu_alloc_GB":peak/2**30,
      "peak_gpu_reserved_GB":torch.cuda.max_memory_reserved()/2**30}
print(json.dumps(summ))
open(a.out+".done","w").write(json.dumps(summ)+"\n")
print("SCORE_DONE")
