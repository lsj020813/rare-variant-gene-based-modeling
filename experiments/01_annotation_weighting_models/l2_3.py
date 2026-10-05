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


import glob, gzip, zipfile, io, json, os, numpy as np
from scipy.stats import spearmanr
W=_config_path('${PROJECT_ROOT}/work'); rng=np.random.default_rng(20260914)
MAP={'tchl':['TC'],'htn':['SBP','DBP'],'dm':['T2D'],'lip':['LDLC','TG','HDLC']}
def swap(k): a=k.split(':'); return f'{a[0]}:{a[1]}:{a[3]}:{a[2]}'
units={}
for T in MAP:
    for f in glob.glob(f'{W}/run_ourfm/our_pip/{T}.chr*.tsv'):
        with open(f) as fh:
            hdr=fh.readline().rstrip('\n').split('\t'); ci={h:i for i,h in enumerate(hdr)}
            for line in fh:
                q=line.rstrip('\n').split('\t')
                if q[ci['cs_id']]=='NA': continue
                units.setdefault((T,q[ci['region']],q[ci['cs_id']]),[]).append((q[ci['variant_hg19']],float(q[ci['pip']])))
units={k:v for k,v in units.items() if len(v)>=2}
need={T:set() for T in MAP}
for (T,r,c),lst in units.items():
    for k,_ in lst: need[T].add(k); need[T].add(swap(k))
bbj={T:{} for T in MAP}; nrows={}
for T,bts in MAP.items():
    for bt in bts:
        zp=f'{W}/ref/bbj_fm/hum0197.v5.finemap.{bt}.v1.zip'
        with zipfile.ZipFile(zp) as z:
            name=[n for n in z.namelist() if 'SuSiE' in n][0]
            with gzip.open(io.BytesIO(z.read(name)),'rt') as fh:
                hdr=fh.readline().rstrip('\n').split('\t'); ic,ip,ia1,ia2,ipip=[hdr.index(x) for x in ('chromosome','position','allele1','allele2','pip')]
                n=0
                for line in fh:
                    f=line.rstrip('\n').split('\t'); k=f"{f[ic].replace('chr','')}:{f[ip]}:{f[ia1]}:{f[ia2]}"
                    n+=1
                    if k in need[T]:
                        try: p=float(f[ipip])
                        except: continue
                        kk=k if k in need[T] and any(k==x for x in (k,)) else k
                        bbj[T][k]=max(p,bbj[T].get(k,-1))
                nrows[bt]=n
def lookup(T,k): 
    v=bbj[T].get(k); 
    if v is None: v=bbj[T].get(swap(k))
    return v
present=0; total=0; rows=[]; rhos=[]; ce_bbj=[]; ce_flat=[]; n_units_used=0
for (T,r,c),lst in sorted(units.items()):
    y=np.array([p for _,p in lst]); b=np.array([lookup(T,k) for k,_ in lst],dtype=float)
    total+=len(lst); present+=int(np.isfinite(b).sum())
    have=np.isfinite(b)
    if have.sum()<2: continue
    yy=y[have]/y[have].sum(); bb=b[have]; n_units_used+=1
    pb=(bb+1e-4)/(bb+1e-4).sum(); ce_bbj.append(float(-(yy*np.log(pb)).sum())); ce_flat.append(float(-(yy*np.log(1.0/len(yy))).sum()))
    if have.sum()>=3 and np.std(bb)>0 and np.std(yy)>0: rhos.append(float(spearmanr(yy,bb).correlation))
obs_d=float(np.mean(ce_flat)-np.mean(ce_bbj)); obs_rho=float(np.median(rhos)) if rhos else None
nd=[]; nr=[]
for it in range(200):
    d=[]; rr=[]
    for (T,r,c),lst in units.items():
        y=np.array([p for _,p in lst]); b=np.array([lookup(T,k) for k,_ in lst],dtype=float); have=np.isfinite(b)
        if have.sum()<2: continue
        yy=y[have]/y[have].sum(); bb=rng.permutation(b[have]); pb=(bb+1e-4)/(bb+1e-4).sum()
        d.append(-(yy*np.log(1.0/len(yy))).sum()+(yy*np.log(pb)).sum())
        if have.sum()>=3 and np.std(bb)>0 and np.std(yy)>0: rr.append(spearmanr(yy,bb).correlation)
    nd.append(np.mean(d)); nr.append(np.median(rr) if rr else np.nan)
res=dict(prereg='v160',bbj_rows=nrows,n_units=len(units),n_units_with_bbj_ge2=n_units_used,n_var=total,frac_var_with_bbj_pip=round(present/total,4),
         ce_flat=float(np.mean(ce_flat)),ce_bbj_pip=float(np.mean(ce_bbj)),dCE_bbjpip_vs_flat=obs_d,perm_p_dCE=float(np.mean(np.array(nd)>=obs_d)),null_q95_dCE=float(np.quantile(nd,.95)),
         n_cs_rho=len(rhos),median_rho=obs_rho,frac_rho_pos=float(np.mean(np.array(rhos)>0)) if rhos else None,perm_p_rho=float(np.nanmean(np.array(nr)>=obs_rho)) if rhos else None,null_q95_rho=float(np.nanquantile(nr,.95)) if rhos else None,
         frac_bbj_pip_gt0p1_among_present=float(np.mean([v>0.1 for T in bbj for v in bbj[T].values()])) if any(bbj.values()) else None)
json.dump(res,open(f'{W}/run_l2/l2_3_result.json','w'),indent=1); print(json.dumps(res)); print('L2_3_DONE')
