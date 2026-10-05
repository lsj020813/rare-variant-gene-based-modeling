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


import sys, os, glob, gzip, json, time, math
sys.path.insert(0, _config_path('${PROJECT_ROOT}/work/run_band15/model_v10_out'))
import l1_train_v10 as L; L.load_libraries(4)
import numpy as np
from scipy.optimize import minimize
from scipy.special import expit
W=_config_path('${PROJECT_ROOT}/work'); TA=f'{W}/run_trackA'; MODEL=f'{W}/run_band15/model_v10_out/e6_runs/primary_f0/phi.json'
model=json.load(open(MODEL)); cols=model['design']['cols']; beta0=np.asarray(model['coefficients']['shared'],float); b0=float(model['default_intercept'])
rng=np.random.default_rng(20260914); SMOKE='--smoke' in sys.argv
def num(v): return float('nan') if v in ('','NA','nan') else float(v)
def swap(k): a=k.split(':'); return f'{a[0]}:{a[1]}:{a[3]}:{a[2]}'
trA={}
with gzip.open(f'{TA}/annot_trA/trA_annot_missing.tsv.gz','rt') as fh:
    hdr=fh.readline().rstrip('\n').split('\t'); ci={c:hdr.index(c) for c in cols}; ki=hdr.index('variant_hg19')
    for line in fh:
        a=line.rstrip('\n').split('\t'); trA[a[ki]]=[num(a[ci[c]]) for c in cols]
bbj={}
for bf in glob.glob(f'{TA}/annot/*.bbj.tsv'):
    with open(bf) as fh:
        hdr=fh.readline().rstrip('\n').split('\t'); ci={c:hdr.index(c) for c in cols}; ki=hdr.index('variant_hg19')
        for line in fh:
            a=line.rstrip('\n').split('\t'); k=a[ki] if a[0]=='0' else swap(a[ki]); bbj[k]=[num(a[ci[c]]) for c in cols]
units={}
nvar=0
for T in ['tchl','htn','dm','lip']:
    for f in glob.glob(f'{W}/run_ourfm/our_pip/{T}.chr*.tsv'):
        with open(f) as fh:
            hdr=fh.readline().rstrip('\n').split('\t'); ci={h:i for i,h in enumerate(hdr)}
            for line in fh:
                q=line.rstrip('\n').split('\t')
                if q[ci['cs_id']]=='NA': continue
                nvar+=1; units.setdefault((T,q[ci['region']],q[ci['cs_id']]),[]).append((q[ci['variant_hg19']],float(q[ci['pip']])))
keys=[];Y=[];cs_idx=[];reg_of=[];src_cnt={'bbj':0,'trA':0,'none':0}; X=[]
u=0
for (T,rid,cs),lst in sorted(units.items()):
    rows=[]
    for k,pip in lst:
        if k in bbj: rows.append((k,pip,bbj[k])); src_cnt['bbj']+=1
        elif k in trA: rows.append((k,pip,trA[k])); src_cnt['trA']+=1
        else: src_cnt['none']+=1
    if len(rows)<2: continue
    for k,pip,x in rows: keys.append(k); Y.append(pip); X.append(x); cs_idx.append(u); reg_of.append(rid)
    u+=1
X=np.asarray(X,float); Y=np.asarray(Y); cs_idx=np.asarray(cs_idx); reg_of=np.asarray(reg_of); NU=u
if 't1_na' in cols:
    na=X[:,cols.index('t1_na')]==1
    for col in L.S1[:7]+['cadd']:
        if col in cols: X[na,cols.index(col)]=np.nan
B=L.transform_design(model['design'],X,cols).astype(np.float64)
Ysum=np.bincount(cs_idx,weights=Y,minlength=NU); Yn=Y/Ysum[cs_idx]
print(json.dumps(dict(n_cs_var=nvar,n_units=NU,n_rows=len(Y),src=src_cnt,n_regions=len(set(reg_of)),n_coef=len(beta0),B_shape=B.shape)),flush=True)
def ce_and_grad(beta,rows,lam,yn):
    eta=B[rows]@beta+b0; phi=expit(eta); c=cs_idx[rows]
    S=np.bincount(c,weights=phi,minlength=NU); p=phi/S[c]
    ncs=len(np.unique(c)); ce=-(yn[rows]*np.log(np.clip(p,1e-12,1))).sum()/ncs
    g_eta=(1-phi)*(p-yn[rows])/ncs
    g=B[rows].T@g_eta+2*lam*(beta-beta0)
    return ce+lam*((beta-beta0)**2).sum(), g
def fit(rows,lam,yn):
    if lam==np.inf: return beta0.copy()
    r=minimize(lambda b: ce_and_grad(b,rows,lam,yn), beta0, jac=True, method='L-BFGS-B', options=dict(maxiter=500))
    return r.x
def ce_only(beta,rows,yn): return ce_and_grad(beta,rows,0.0,yn)[0]
def ce_flat(rows,yn):
    c=cs_idx[rows]; M=np.bincount(c,minlength=NU)[c]; ncs=len(np.unique(c)); return -(yn[rows]*np.log(1.0/M)).sum()/ncs
LAMS=[np.inf,100,30,10,3,1,0.3,0.1,0.0]
regions=np.array(sorted(set(reg_of))); folds=np.array_split(rng.permutation(regions),5)
def nested(yn, do_wise=True):
    out=[]
    for k in range(5):
        te=np.where(np.isin(reg_of,folds[k]))[0]; tr=np.where(~np.isin(reg_of,folds[k]))[0]
        if len(te)==0: continue
        treg=np.array(sorted(set(reg_of[tr]))); ifolds=np.array_split(rng.permutation(treg),4); inner={lam:0.0 for lam in LAMS}
        for j in range(4):
            ite=tr[np.isin(reg_of[tr],ifolds[j])]; itr=tr[~np.isin(reg_of[tr],ifolds[j])]
            if len(ite)==0: continue
            for lam in LAMS: inner[lam]+=ce_only(fit(itr,lam,yn),ite,yn)
        lam_hat=min(LAMS,key=lambda l: inner[l])
        curve={str(lam): float(ce_only(fit(tr,lam,yn),te,yn)) for lam in LAMS}
        rec=dict(fold=k,n_test_cs=int(len(np.unique(cs_idx[te]))),lam_hat=(None if lam_hat==np.inf else lam_hat),ce_flat=float(ce_flat(te,yn)),ce_curve=curve,ce_hat=curve[str(lam_hat)])
        if do_wise:
            bfull=fit(tr,0.0,yn); rec['wise']={str(a): float(ce_only(a*bfull+(1-a)*beta0,te,yn)) for a in np.round(np.arange(0,1.01,0.1),1)}
            d=bfull-beta0; names=model['design']['coef_names']; top=np.argsort(-np.abs(d))[:5]; rec['top_shift']=[(names[i],round(float(d[i]),4)) for i in top]
        out.append(rec)
    return out
t0=time.time(); obs=nested(Yn); print('obs done',round(time.time()-t0),flush=True)
def agg(res):
    f=np.mean([r['ce_flat'] for r in res]); inf=np.mean([r['ce_curve']['inf'] for r in res]); hat=np.mean([r['ce_hat'] for r in res])
    return dict(ce_flat=float(f),ce_bbj=float(inf),ce_hat=float(hat),dCE_transfer=float(f-inf),dCE_finetune=float(inf-hat))
A=agg(obs)
NP=20 if SMOKE else 200; null_tr=[]; null_ft=[]
for b in range(NP):
    yp=Yn.copy()
    for uu in range(NU):
        idx=np.where(cs_idx==uu)[0]; yp[idx]=rng.permutation(yp[idx])
    r=nested(yp,do_wise=False); a=agg(r); null_tr.append(a['dCE_transfer']); null_ft.append(a['dCE_finetune'])
    if b%20==19: print('perm',b+1,round(time.time()-t0),flush=True)
res=dict(prereg='v158',n_units=NU,n_rows=len(Y),src=src_cnt,folds=obs,agg=A,
         curve_mean={str(l):float(np.mean([r['ce_curve'][str(l)] for r in obs])) for l in LAMS},
         wise_mean={a:float(np.mean([r['wise'][a] for r in obs])) for a in obs[0]['wise']},
         lam_hat_by_fold=[r['lam_hat'] for r in obs],
         perm=dict(n=NP,transfer_p=float(np.mean(np.array(null_tr)>=A['dCE_transfer'])),transfer_q95=float(np.quantile(null_tr,.95)),
                   finetune_p=float(np.mean(np.array(null_ft)>=A['dCE_finetune'])),finetune_q95=float(np.quantile(null_ft,.95))),
         elapsed_s=round(time.time()-t0))
os.makedirs(f'{W}/run_l2',exist_ok=True); json.dump(res,open(f'{W}/run_l2/l2_2_result{"_smoke" if SMOKE else ""}.json','w'),indent=1)
print(json.dumps({k:v for k,v in res.items() if k!='folds'},default=str)); print('L2_2_DONE')
