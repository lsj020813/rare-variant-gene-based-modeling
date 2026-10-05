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


import sys, os, glob, gzip, json, time
sys.path.insert(0, _config_path('${PROJECT_ROOT}/work/run_band15/model_v10_out'))
import l1_train_v10 as L; L.load_libraries(4)
import numpy as np
from scipy.optimize import minimize
from scipy.special import expit
W=_config_path('${PROJECT_ROOT}/work'); TA=f'{W}/run_trackA'; BL=f'{W}/ref/annot/bbj_l1'
MODEL=f'{W}/run_band15/model_v10_out/e6_runs/primary_f0/phi.json'
model=json.load(open(MODEL)); cols=model['design']['cols']; beta0=np.asarray(model['coefficients']['shared'],float); b0=float(model['default_intercept'])
names=model['design']['coef_names']; rng=np.random.default_rng(20260914); SMOKE='--smoke' in sys.argv
BBJ_TRAITS=['TC','SBP','DBP','T2D','LDLC','TG','HDLC']
def num(v): return float('nan') if v in ('','NA','nan') else float(v)
def swap(k): a=k.split(':'); return f'{a[0]}:{a[1]}:{a[3]}:{a[2]}'
bunits={}
for T in BBJ_TRAITS:
    with gzip.open(f'{BL}/bbj_pip/{T}.tsv.gz','rt') as fh:
        hdr=fh.readline().rstrip('\n').split('\t'); ci={h:i for i,h in enumerate(hdr)}
        for line in fh:
            q=line.rstrip('\n').split('\t'); cid=int(q[ci['cs_id']])
            if cid<0: continue
            bunits.setdefault((T,q[ci['region']],cid),[]).append((q[ci['variant_hg19']],float(q[ci['pip']])))
bunits={k:v for k,v in bunits.items() if len(v)>=2}
need=set(k for lst in bunits.values() for k,_ in lst)
bann={}
with gzip.open(f'{BL}/bbj_annot.tsv.gz','rt') as fh:
    hdr=fh.readline().rstrip('\n').split('\t'); ci={c:hdr.index(c) for c in cols}; ki=hdr.index('variant_hg19')
    for line in fh:
        a=line.split('\t'); k=a[ki]
        if k in need: bann[k]=[num(a[ci[c]].strip()) for c in cols]
trA={}
with gzip.open(f'{TA}/annot_trA/trA_annot_missing.tsv.gz','rt') as fh:
    hdr=fh.readline().rstrip('\n').split('\t'); ci={c:hdr.index(c) for c in cols}; ki=hdr.index('variant_hg19')
    for line in fh:
        a=line.rstrip('\n').split('\t'); trA[a[ki]]=[num(a[ci[c]]) for c in cols]
kann={}
for bf in glob.glob(f'{TA}/annot/*.bbj.tsv'):
    with open(bf) as fh:
        hdr=fh.readline().rstrip('\n').split('\t'); ci={c:hdr.index(c) for c in cols}; ki=hdr.index('variant_hg19')
        for line in fh:
            a=line.rstrip('\n').split('\t'); k=a[ki] if a[0]=='0' else swap(a[ki]); kann[k]=[num(a[ci[c]]) for c in cols]
kunits={}
for T in ['tchl','htn','dm','lip']:
    for f in glob.glob(f'{W}/run_ourfm/our_pip/{T}.chr*.tsv'):
        with open(f) as fh:
            hdr=fh.readline().rstrip('\n').split('\t'); ci={h:i for i,h in enumerate(hdr)}
            for line in fh:
                q=line.rstrip('\n').split('\t')
                if q[ci['cs_id']]=='NA': continue
                kunits.setdefault((T,q[ci['region']],q[ci['cs_id']]),[]).append((q[ci['variant_hg19']],float(q[ci['pip']])))
X=[];Y=[];cs=[];dom=[];reg=[];u=0
def add(units,ann,d):
    global u
    for key,lst in sorted(units.items()):
        rows=[(k,p,ann[k]) for k,p in lst if k in ann]
        if len(rows)<2: continue
        for k,p,x in rows: X.append(x);Y.append(p);cs.append(u);dom.append(d);reg.append(f'{d}:{key[1]}')
        u+=1
add(bunits,bann,0); nB=u; add(kunits,kann,1); NU=u
X=np.asarray(X,float);Y=np.asarray(Y);cs=np.asarray(cs);dom=np.asarray(dom);reg=np.asarray(reg)
if 't1_na' in cols:
    na=X[:,cols.index('t1_na')]==1
    for col in L.S1[:7]+['cadd']:
        if col in cols: X[na,cols.index(col)]=np.nan
B=L.transform_design(model['design'],X,cols).astype(np.float64); q=B.shape[1]
Z=np.zeros((len(Y),3*q)); Z[:,:q]=B; Z[dom==0,q:2*q]=B[dom==0]; Z[dom==1,2*q:]=B[dom==1]
Ysum=np.bincount(cs,weights=Y,minlength=NU); Yn=Y/Ysum[cs]
nK=NU-nB; wK=nB/max(nK,1)
print(json.dumps(dict(bbj_cs=nB,kor_cs=nK,rows=int(len(Y)),rows_bbj=int((dom==0).sum()),rows_kor=int((dom==1).sum()),q=q,w_kor=round(wK,2))),flush=True)
theta0=np.concatenate([beta0,np.zeros(q),np.zeros(q)])
def loss(theta,rows,lam,yn,w):
    eta=Z[rows]@theta+b0; phi=expit(eta); c=cs[rows]
    S=np.bincount(c,weights=phi,minlength=NU); p=phi/S[c]
    wr=np.where(dom[rows]==1,w,1.0); ncs_eff=wr[np.unique(c,return_index=True)[1]].sum()
    ce=-(wr*yn[rows]*np.log(np.clip(p,1e-12,1))).sum()/ncs_eff
    g_eta=wr*(1-phi)*(p-yn[rows])/ncs_eff
    bc,bs,bk=theta[:q],theta[q:2*q],theta[2*q:]
    pen=lam*(((bc-beta0)**2).sum()+(bs**2).sum()+(bk**2).sum())
    g=Z[rows].T@g_eta+2*lam*np.concatenate([bc-beta0,bs,bk])
    return ce+pen,g
def fit(rows,lam,yn,w):
    return minimize(lambda t: loss(t,rows,lam,yn,w), theta0, jac=True, method='L-BFGS-B', options=dict(maxiter=600)).x
def ce_dom(theta,rows,yn,d):
    r=rows[dom[rows]==d]
    if len(r)==0: return None
    eta=Z[r]@theta+b0; phi=expit(eta); c=cs[r]; S=np.bincount(c,weights=phi,minlength=NU); p=phi/S[c]
    return float(-(yn[r]*np.log(np.clip(p,1e-12,1))).sum()/len(np.unique(c)))
LAMS=[10,3,1,0.3,0.1]
regions=np.array(sorted(set(reg))); folds=np.array_split(rng.permutation(regions),5)
def run(yn,w,lam_fixed=None):
    out=[]
    for k in range(5):
        te=np.where(np.isin(reg,folds[k]))[0]; tr=np.where(~np.isin(reg,folds[k]))[0]
        if lam_fixed is None:
            treg=np.array(sorted(set(reg[tr]))); ifolds=np.array_split(rng.permutation(treg),3); sc={l:0.0 for l in LAMS}
            for j in range(3):
                ite=tr[np.isin(reg[tr],ifolds[j])]; itr=tr[~np.isin(reg[tr],ifolds[j])]
                for l in LAMS:
                    th=fit(itr,l,yn,w); a=ce_dom(th,ite,yn,0); b=ce_dom(th,ite,yn,1); sc[l]+=(a or 0)+(b or 0)
            lam=min(LAMS,key=lambda l: sc[l])
        else: lam=lam_fixed
        th=fit(tr,lam,yn,w)
        out.append(dict(fold=k,lam=lam,theta=th,ce_bbj=ce_dom(th,te,yn,0),ce_kor=ce_dom(th,te,yn,1)))
    return out
t0=time.time(); obs=run(Yn,wK); print('obs done',round(time.time()-t0),flush=True)
TH=np.array([o['theta'] for o in obs]); bk=TH[:,2*q:]; bs=TH[:,q:2*q]
def stable(M): return (np.abs(np.sign(M).sum(0))>=4)
mk=np.abs(bk.mean(0)); ms=np.abs(bs.mean(0)); lam_hat=[o['lam'] for o in obs]
NP=10 if SMOKE else 200; null_k=[]; null_s=[]
lam_mode=max(set(lam_hat),key=lam_hat.count)
for b in range(NP):
    yp=Yn.copy()
    for uu in range(NU):
        idx=np.where(cs==uu)[0]; yp[idx]=rng.permutation(yp[idx])
    r=run(yp,wK,lam_fixed=lam_mode); T2=np.array([o['theta'] for o in r])
    null_k.append(np.abs(T2[:,2*q:].mean(0))); null_s.append(np.abs(T2[:,q:2*q].mean(0)))
    if b%20==19: print('perm',b+1,round(time.time()-t0),flush=True)
null_k=np.array(null_k); null_s=np.array(null_s)
qk=np.quantile(null_k,.95,axis=0); qs=np.quantile(null_s,.95,axis=0)
sigK=[names[i] for i in range(q) if stable(bk)[i] and mk[i]>qk[i]]; sigS=[names[i] for i in range(q) if stable(bs)[i] and ms[i]>qs[i]]
res=dict(prereg='v162',bbj_cs=nB,kor_cs=nK,rows=int(len(Y)),w_kor=round(wK,2),lam_hat_by_fold=lam_hat,
         ce_heldout=dict(bbj=float(np.mean([o['ce_bbj'] for o in obs])),kor=float(np.mean([o['ce_kor'] for o in obs if o['ce_kor'] is not None]))),
         k_KOR=len(sigK),sig_KOR=sigK,k_BBJ=len(sigS),sig_BBJ=sigS,
         top_KOR_by_mean=[(names[i],round(float(bk.mean(0)[i]),4),round(float(qk[i]),4)) for i in np.argsort(-mk)[:8]],
         perm_n=NP,elapsed_s=round(time.time()-t0))
json.dump(res,open(f'{W}/run_l2/l2_4_result{"_smoke" if SMOKE else ""}.json','w'),indent=1); print(json.dumps(res)); print('L2_4_DONE')
