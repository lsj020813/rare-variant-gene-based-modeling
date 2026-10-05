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


import sys, os, glob, gzip, zipfile, io, json, time
sys.path.insert(0, _config_path('${PROJECT_ROOT}/work/run_band15/model_v10_out'))
import l1_train_v10 as L; L.load_libraries(4)
import numpy as np
from scipy.optimize import minimize
from scipy.special import expit
from scipy.stats import spearmanr
W=_config_path('${PROJECT_ROOT}/work'); K=f'{W}/run_kcps2fm'; TA=f'{W}/run_trackA'
MODEL=f'{W}/run_band15/model_v10_out/e6_runs/primary_f0/phi.json'
model=json.load(open(MODEL)); cols=model['design']['cols']; beta0=np.asarray(model['coefficients']['shared'],float); b0=float(model['default_intercept'])
names=model['design']['coef_names']; rng=np.random.default_rng(20260915); SMOKE='--smoke' in sys.argv; INTERIM='--interim' in sys.argv; NOPERM='--noperm' in sys.argv
PURITY=0.5; LAM_MAX=0.05
BBJMAP={'tchl':['TC'],'htn':['SBP','DBP'],'dm':['T2D'],'lip':['LDLC','TG','HDLC']}
if INTERIM: BBJMAP={k:v for k,v in BBJMAP.items() if k!='lip'}
def num(v): return float('nan') if v in ('','NA','nan') else float(v)
def swap(k): a=k.split(':'); return f'{a[0]}:{a[1]}:{a[3]}:{a[2]}'
units={}; stats=dict(regions=0, cs_raw=0, cs_purity_drop=0, regions_lambda_drop=0, cs_size1=0)
for T in BBJMAP:
    for pf in glob.glob(f'{K}/fm/{T}/*.pip.tsv'):
        rid=os.path.basename(pf)[:-8]
        if not os.path.exists(f'{K}/fm/{T}/{rid}.done'): continue
        stats['regions']+=1
        lam=None
        for line in open(f'{K}/fm/{T}/{rid}.susie.log'):
            if line.startswith('SUSIE_OK'):
                for tok in line.split():
                    if tok.startswith('lambda='): lam=float(tok[7:])
        if lam is not None and lam>LAM_MAX: stats['regions_lambda_drop']+=1; continue
        pur={}
        sf=f'{K}/fm/{T}/{rid}.summary.tsv'
        if os.path.exists(sf):
            with open(sf) as fh:
                hdr=fh.readline().rstrip('\n').split('\t'); ci={h:i for i,h in enumerate(hdr)}; q=fh.readline().rstrip('\n').split('\t')
                if 'cs_min_abs_corr' in ci and len(q)>ci['cs_min_abs_corr']:
                    import re as _re
                    vals=[float(x) for x in _re.split(r'[;,|]',q[ci['cs_min_abs_corr']]) if x.strip() not in ('','NA')]
                    for j,v in enumerate(vals): pur[str(j+1)]=v; pur[f'L{j+1}']=v
        with open(pf) as fh:
            hdr=fh.readline().rstrip('\n').split('\t'); ci={h:i for i,h in enumerate(hdr)}
            for line in fh:
                q=line.rstrip('\n').split('\t'); c=q[ci['cs_id']]
                if c in ('NA','-1',''): continue
                units.setdefault((T,rid,c),[]).append((q[ci['variant_hg19']],float(q[ci['pip']])))
        for c in set(k[2] for k in units if k[1]==rid):
            stats['cs_raw']+=1
            if c in pur and pur[c]<PURITY: stats['cs_purity_drop']+=1; units.pop((T,rid,c),None)
units={k:v for k,v in units.items() if len(v)>=2}; stats['cs_size1']=stats['cs_raw']-stats['cs_purity_drop']-len(units)
ann={}
with gzip.open(f'{TA}/annot_trA/trA_annot_missing.tsv.gz','rt') as fh:
    hdr=fh.readline().rstrip('\n').split('\t'); ci={c:hdr.index(c) for c in cols}; ki=hdr.index('variant_hg19')
    for line in fh:
        a=line.rstrip('\n').split('\t'); ann[a[ki]]=[num(a[ci[c]]) for c in cols]
for bf in glob.glob(f'{TA}/annot/*.bbj.tsv'):
    with open(bf) as fh:
        hdr=fh.readline().rstrip('\n').split('\t'); ci={c:hdr.index(c) for c in cols}; ki=hdr.index('variant_hg19')
        for line in fh:
            a=line.rstrip('\n').split('\t'); k=a[ki] if a[0]=='0' else swap(a[ki]); ann.setdefault(k,[num(a[ci[c]]) for c in cols])
need=set(k for lst in units.values() for k,_ in lst)-set(ann)
if need:
    with gzip.open(f'{W}/ref/annot/bbj_l1/bbj_annot.tsv.gz','rt') as fh:
        hdr=fh.readline().rstrip('\n').split('\t'); have=[c for c in cols if c in hdr]; ci={c:hdr.index(c) for c in have}; ki=hdr.index('variant_hg19')
        for line in fh:
            a=line.split('\t'); k=a[ki]
            if k in need or swap(k) in need:
                row=[num(a[ci[c]].strip()) if c in ci else float('nan') for c in cols]; ann[k]=row; ann[swap(k)]=row
bbjpip={T:{} for T in BBJMAP}; needk={T:set() for T in BBJMAP}
for (T,r,c),lst in units.items():
    for k,_ in lst: needk[T].add(k); needk[T].add(swap(k))
for T,bts in BBJMAP.items():
    for bt in bts:
        zp=f'{W}/ref/bbj_fm/hum0197.v5.finemap.{bt}.v1.zip'
        with zipfile.ZipFile(zp) as z:
            name=[n for n in z.namelist() if 'SuSiE' in n][0]
            with gzip.open(io.BytesIO(z.read(name)),'rt') as fh:
                hdr=fh.readline().rstrip('\n').split('\t'); ic,ip,ia1,ia2,ipip=[hdr.index(x) for x in ('chromosome','position','allele1','allele2','pip')]
                for line in fh:
                    f=line.rstrip('\n').split('\t'); k=f"{f[ic].replace('chr','')}:{f[ip]}:{f[ia1]}:{f[ia2]}"
                    if k in needk[T]:
                        try: bbjpip[T][k]=max(float(f[ipip]),bbjpip[T].get(k,-1))
                        except: pass
def bbjlook(T,k):
    v=bbjpip[T].get(k); return v if v is not None else bbjpip[T].get(swap(k))
X=[];Y=[];cs=[];reg=[];trait=[];BP=[];u=0; src=dict(annot=0,noannot=0)
for (T,rid,c),lst in sorted(units.items()):
    rows=[(k,p,ann[k]) for k,p in lst if k in ann]; src['annot']+=len(rows); src['noannot']+=len(lst)-len(rows)
    if len(rows)<2: continue
    for k,p,x in rows: X.append(x);Y.append(p);cs.append(u);reg.append(rid);trait.append(T);BP.append(bbjlook(T,k))
    u+=1
X=np.asarray(X,float);Y=np.asarray(Y);cs=np.asarray(cs);reg=np.asarray(reg);trait=np.asarray(trait);NU=u
BPa=np.array([np.nan if v is None else v for v in BP],float)
if 't1_na' in cols:
    na=X[:,cols.index('t1_na')]==1
    for col in L.S1[:7]+['cadd']:
        if col in cols: X[na,cols.index(col)]=np.nan
B=L.transform_design(model['design'],X,cols).astype(np.float64)
Ysum=np.bincount(cs,weights=Y,minlength=NU); Yn=Y/Ysum[cs]
print(json.dumps(dict(stats=stats,n_units=NU,n_rows=int(len(Y)),src=src,by_trait={T:int(len(set(cs[trait==T]))) for T in BBJMAP},bbj_pip_coverage=float(np.isfinite(BPa).mean()))),flush=True)
def cs_ce(p,yn,rows):
    c=cs[rows]; S=np.bincount(c,weights=p,minlength=NU); q=np.clip(p/S[c],1e-12,1)
    per=np.bincount(c,weights=-yn[rows]*np.log(q),minlength=NU); return per
def cs_flat(yn,rows):
    c=cs[rows]; M=np.bincount(c,minlength=NU); return np.bincount(c,weights=yn[rows]*np.log(M[c]),minlength=NU)
def cs_H(yn,rows):
    c=cs[rows]; yy=np.clip(yn[rows],1e-12,1); return np.bincount(c,weights=-yn[rows]*np.log(yy),minlength=NU)
def phi_of(beta,rows): return expit(B[rows]@beta+b0)
def agg(per_flat,per_phi,per_H,units_idx):
    if len(units_idx)==0: return dict(dCE=None,SS=None,n_cs=0,n_cs_informative=0,den_sum=0.0)
    d=per_flat[units_idx]-per_phi[units_idx]; den=per_flat[units_idx]-per_H[units_idx]
    ok=den>1e-9
    return dict(dCE=float(d.mean()), SS=float((d[ok]).sum()/den[ok].sum()) if ok.any() else None, n_cs=int(len(units_idx)), n_cs_informative=int(ok.sum()), den_sum=float(den[ok].sum()))
allrows=np.arange(len(Y)); allunits=np.arange(NU)
def U(T): return np.array(sorted(set(cs[trait==T])),dtype=int)
pf=cs_flat(Yn,allrows); pH=cs_H(Yn,allrows)
res=dict(prereg='v164/v168',stats=stats,n_units=NU,n_rows=int(len(Y)),src=src)
pphi=cs_ce(phi_of(beta0,allrows),Yn,allrows)
res['transfer']={T if T else 'ALL': agg(pf,pphi,pH, allunits if not T else U(T)) for T in ['']+list(BBJMAP)}
hb=np.isfinite(BPa); rows_b=allrows[hb]; units_b=np.array(sorted(set(cs[rows_b])),dtype=int)
if len(units_b):
    pb=BPa.copy(); pb[~hb]=0; pb=pb+1e-4; pBB=cs_ce(pb,Yn,allrows)
    res['bbj_pip_as_pred']=agg(pf,pBB,pH,units_b); res['bbj_pip_as_pred']['pip_gt0p1_frac']=float((BPa[hb]>0.1).mean())
    rhos=[]
    for uu in units_b:
        idx=allrows[(cs==uu)&hb]
        if len(idx)>=3 and np.std(BPa[idx])>0 and np.std(Yn[idx])>0: rhos.append(spearmanr(Yn[idx],BPa[idx]).correlation)
    res['rho_bbj_vs_kcps2']=dict(n_cs=len(rhos),median=float(np.median(rhos)) if rhos else None,frac_pos=float(np.mean(np.array(rhos)>0)) if rhos else None)
res['label']=dict(pip_gt0p1_frac=float((Y>0.1).mean()),pip_gt0p5_frac=float((Y>0.5).mean()),median_cs_size=int(np.median(np.bincount(cs))), H_over_logM_median=float(np.median(pH/np.log(np.bincount(cs)))))
def loss(beta,rows,lam):
    eta=B[rows]@beta+b0; phi=expit(eta); c=cs[rows]; S=np.bincount(c,weights=phi,minlength=NU); p=phi/S[c]
    ncs=len(np.unique(c)); ce=-(Yn[rows]*np.log(np.clip(p,1e-12,1))).sum()/ncs
    g=B[rows].T@((1-phi)*(p-Yn[rows])/ncs)+2*lam*(beta-beta0)
    return ce+lam*((beta-beta0)**2).sum(), g
def fit(rows,lam):
    if lam==np.inf: return beta0.copy()
    return minimize(lambda b: loss(b,rows,lam), beta0, jac=True, method='L-BFGS-B', options=dict(maxiter=500)).x
LAMS=[np.inf,100,30,10,3,1,0.3,0.1,0.0]
regions=np.array(sorted(set(reg))); folds=np.array_split(rng.permutation(regions),5)
def nested(yn):
    out=[]
    for k in range(5):
        te=np.where(np.isin(reg,folds[k]))[0]; tr=np.where(~np.isin(reg,folds[k]))[0]
        if len(te)==0: continue
        treg=np.array(sorted(set(reg[tr]))); ifolds=np.array_split(rng.permutation(treg),4); inner={l:0.0 for l in LAMS}
        for j in range(4):
            ite=tr[np.isin(reg[tr],ifolds[j])]; itr=tr[~np.isin(reg[tr],ifolds[j])]
            if len(ite)==0: continue
            for l in LAMS: inner[l]+=loss(fit(itr,l),ite,0.0)[0]
        lh=min(LAMS,key=lambda l: inner[l]); ute=np.array(sorted(set(cs[te])),dtype=int)
        bh=fit(tr,lh); pph=cs_ce(phi_of(bh,allrows),yn,allrows); pp0=cs_ce(phi_of(beta0,allrows),yn,allrows)
        out.append(dict(fold=k,lam_hat=(None if lh==np.inf else lh),n_test_cs=int(len(ute)),
                        ce_flat=float(pf[ute].mean()),ce_bbj=float(pp0[ute].mean()),ce_hat=float(pph[ute].mean()),
                        SS_bbj=agg(pf,pp0,pH,ute)['SS'],SS_hat=agg(pf,pph,pH,ute)['SS']))
    return out
t0=time.time(); ft=nested(Yn); res['finetune']=dict(folds=ft,lam_hat=[f['lam_hat'] for f in ft],
    ce_flat=float(np.mean([f['ce_flat'] for f in ft])),ce_bbj=float(np.mean([f['ce_bbj'] for f in ft])),ce_hat=float(np.mean([f['ce_hat'] for f in ft])),
    SS_bbj=float(np.mean([f['SS_bbj'] for f in ft if f['SS_bbj'] is not None])),SS_hat=float(np.mean([f['SS_hat'] for f in ft if f['SS_hat'] is not None])))
print('obs done',round(time.time()-t0),flush=True)
NP=0 if NOPERM else (20 if SMOKE else 200); nd=[];ns=[];ndf=[];nsf=[]
phi0=phi_of(beta0,allrows)
for b in range(NP):
    pp=phi0.copy()
    for uu in range(NU):
        idx=np.where(cs==uu)[0]; pp[idx]=rng.permutation(pp[idx])
    a=agg(pf,cs_ce(pp,Yn,allrows),pH,allunits); nd.append(a['dCE']); ns.append(a['SS'])
    if b%50==49: print('perm',b+1,round(time.time()-t0),flush=True)
ns=[x for x in ns if x is not None]
if NP==0: res['perm']=None
else: res['perm']=dict(n=NP,dCE_p=float(np.mean(np.array(nd)>=res['transfer']['ALL']['dCE'])),dCE_q95=float(np.quantile(nd,.95)),
                 SS_p=float(np.mean(np.array(ns)>=res['transfer']['ALL']['SS'])) if ns else None,SS_q95=float(np.quantile(ns,.95)) if ns else None)
res['elapsed_s']=round(time.time()-t0)
os.makedirs(f'{W}/run_l2',exist_ok=True); res['mode']=('smoke' if SMOKE else ('interim_tchl_htn_dm' if INTERIM else 'full'))+('_noperm' if NOPERM else ''); json.dump(res,open(f'{W}/run_l2/l2_5_result{"_smoke" if SMOKE else ("_interim" if INTERIM else "")}{"_noperm" if NOPERM else ""}.json','w'),indent=1,default=str)
print(json.dumps({k:v for k,v in res.items() if k not in ('finetune',)},default=str)); print('FT',json.dumps({k:v for k,v in res['finetune'].items() if k!='folds'})); print('L2_5_DONE')
