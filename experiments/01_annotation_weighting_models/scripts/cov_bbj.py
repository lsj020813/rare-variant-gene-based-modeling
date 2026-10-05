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


import glob, gzip, os, re, json, sys, time
W=_config_path('${PROJECT_ROOT}/work/run_ourfm'); A=_config_path('${PROJECT_ROOT}/work/ref/annot/bbj_l1/bbj_annot.tsv.gz')
OUT=_config_path('${PROJECT_ROOT}/work/run_trackA'); os.makedirs(OUT+'/annot', exist_ok=True)
t0=time.time()
regions={}
for t in ['dm','htn','lip','tchl']:
    for f in sorted(glob.glob(f'{W}/fm/{t}/*.ld.vars')):
        rid=os.path.basename(f)[:-len('.ld.vars')]
        if '.mi1000' in rid: continue
        v=[l.strip() for l in open(f) if l.strip()]
        regions[(t,rid)]=v
print('n_regions',len(regions),'n_var_total',sum(len(v) for v in regions.values()))
pat=re.compile(r'^\d+:\d+:[ACGT]+:[ACGT]+$')
bad=sum(1 for v in regions.values() for k in v if not pat.match(k)); print('key_format_bad',bad)
want={}
def swap(k):
    a=k.split(':'); return f'{a[0]}:{a[1]}:{a[3]}:{a[2]}'
for tr,v in regions.items():
    for k in v:
        want.setdefault(k,[]).append((tr,0)); want.setdefault(swap(k),[]).append((tr,1))
print('n_unique_keys',len(set(k for v in regions.values() for k in v)))
hit={tr:[0,0] for tr in regions}; n=0; hdr=None
outs={}
with gzip.open(A,'rt') as fh:
    hdr=fh.readline().rstrip('\n').split('\t'); ki=hdr.index('variant_hg19')
    for line in fh:
        n+=1
        k=line.split('\t',1)[0]
        if k in want:
            for tr,sw in want[k]:
                hit[tr][sw]+=1
                if tr not in outs: outs[tr]=open(f'{OUT}/annot/{tr[1]}.bbj.tsv','w'); outs[tr].write('\t'.join(['orient']+hdr)+'\n')
                outs[tr].write(f'{sw}\t'+line)
for o in outs.values(): o.close()
print('bbj_rows',n,'elapsed',round(time.time()-t0))
rows=[]
for tr,v in regions.items():
    rows.append(dict(trait=tr[0],region=tr[1],n_fm=len(v),hit_direct=hit[tr][0],hit_swap=hit[tr][1],cov=round((hit[tr][0]+hit[tr][1])/len(v),4)))
json.dump(dict(header=hdr,rows=rows,bbj_rows=n,key_format_bad=bad),open(f'{OUT}/coverage.json','w'),indent=1)
tot=sum(r['n_fm'] for r in rows); hd=sum(r['hit_direct'] for r in rows); hs=sum(r['hit_swap'] for r in rows)
print('TOTAL n_fm',tot,'hit_direct',hd,'hit_swap',hs,'cov',round((hd+hs)/tot,4))
print('cov_min',min(r['cov'] for r in rows),'cov_median',sorted(r['cov'] for r in rows)[len(rows)//2],'cov_max',max(r['cov'] for r in rows))
open(f'{OUT}/coverage.done','w').write('ok\n')
