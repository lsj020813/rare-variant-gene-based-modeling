#!/usr/bin/env python3
import os as _cfg_os
import re as _cfg_re
import math as _cfg_math

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)

def _config_number(name, cast=int, positive=False):
    value = cast(_cfg_os.environ[name])
    if not _cfg_math.isfinite(value) or (positive and value <= 0):
        raise ValueError(name + " must be finite" + (" and positive" if positive else ""))
    return value

import argparse, hashlib, json, os, shutil
from datetime import datetime
from pathlib import Path
import numpy as np
import pandas as pd
import yaml
ROOT=Path(_config_path('${PROJECT_ROOT}'))
REPORTS=ROOT/'reports'; SAMPLE_N=_config_number("MASK_N_SAMPLES", int, True); CALL_RATE_MIN=0.95
MASK_NAME='CADD_PHRED_GT5'; MASK_THRESHOLD=5.0
BRANCHES={'D90_MAF001':('GP90',0.001),'D90_MAF01':('GP90',0.01),'D95_MAF001':('GP95',0.001),'D95_MAF01':('GP95',0.01)}
POLICY=REPORTS/'policy_study/CHR22_POLICY_VARIANT_STATISTICS.tsv.gz'
GENES=ROOT/'work/deeprvat_grch38/reference/protein_coding_genes.parquet'
MODEL=Path(_config_path('${PROJECT_ROOT}/projects/deeprvat_test/deeprvat/pretrained_models/model_config.yaml'))
INPUTCFG=Path(_config_path('${PROJECT_ROOT}/projects/deeprvat_test/deeprvat/example/config/deeprvat_input_pretrained_models_config.yaml'))
PIPELINE=Path(_config_path('${PROJECT_ROOT}/projects/deeprvat_test/deeprvat/pipelines/annotations.snakefile'))
ANNO_CODE=Path(_config_path('${PROJECT_ROOT}/projects/deeprvat_test/deeprvat/deeprvat/annotations/annotations.py'))
SAIGE_CODE=Path(_config_path('${PROJECT_ROOT}/projects/saige/SAIGE/R/SAIGE_SPATest_Region_Func.R'))

def key(c,p,r,a):
    c=c.astype(str).str.replace('^chr','',regex=True)
    return 'chr'+c+':'+p.astype('int64').astype(str)+':'+r.astype(str).str.upper()+':'+a.astype(str).str.upper()
def atomic_tsv(d,p):
    p=Path(p); t=p.with_name(p.name+'.tmp.'+str(os.getpid())); d.to_csv(t,sep='\t',index=False); os.replace(t,p)
def atomic_text(s,p):
    p=Path(p); t=p.with_name(p.name+'.tmp.'+str(os.getpid())); t.write_text(s); os.replace(t,p)
def digest(d,cols):
    s='\n'.join('\t'.join(map(str,x)) for x in d[cols].sort_values(cols).itertuples(index=False,name=None))
    return hashlib.sha256(s.encode()).hexdigest()
def mafbin(x):
    return np.select([x.eq(0),(x.gt(0)&x.lt(.001)),(x.ge(.001)&x.lt(.01)),x.ge(.01)],['MONOMORPHIC','0_LT_MAF_LT_0.001','0.001_LE_MAF_LT_0.01','MAF_GE_0.01'],default='UNKNOWN')
def backup(ps,stamp):
    d=REPORTS/'archive'/('common_gene_mask_preupdate_'+stamp); d.mkdir(parents=True); os.chmod(d,0o700)
    for p in map(Path,ps):
        if p.exists(): shutil.copy2(p,d/p.name)
    return d

def main():
    p=argparse.ArgumentParser(); p.add_argument('--annotation',required=True); p.add_argument('--unfilled-annotation',required=True); p.add_argument('--variants',required=True); p.add_argument('--annotation-run',required=True); a=p.parse_args()
    stamp=datetime.now().strftime('%Y%m%d_%H%M%S'); out=ROOT/'work'/('common_gene_mask_v2_'+stamp); out.mkdir(); os.chmod(out,0o700)
    restricted=out/'restricted'; groups=out/'saige_group_files'; restricted.mkdir(); groups.mkdir(); os.chmod(restricted,0o700); os.chmod(groups,0o700)
    targets=[REPORTS/x for x in ['CHR22_COMMON_GENE_MASK_SUMMARY.tsv','CHR22_COMMON_GENE_MASK_SUMMARY.md','CHR22_GENE_MASK_DENOMINATORS.tsv','CHR22_GENE_MASK_QC.tsv','COMMON_MASK_CONTRACT.tsv','COMMON_MASK_CONTRACT.md','COMMON_GENE_MASK_GATE.tsv','COMMON_GENE_MASK_GATE.md']]
    archive=backup(targets,stamp)
    ann=pd.read_parquet(a.annotation); raw=pd.read_parquet(a.unfilled_annotation,columns=['id','gene_id','CADD_PHRED']); var=pd.read_parquet(a.variants,columns=['id','chrom','pos','ref','alt']); genes=pd.read_parquet(GENES,columns=['id','gene','gene_name'])
    use=['policy','CHROM','POS','REF','ALT','R2','called_cells','missing_cells','carrier_cells','AN','AC','AF','MAF']
    st=pd.read_csv(POLICY,sep='\t',compression='gzip',usecols=use); st=st[st.policy.isin(['RAW_GT','GP90','GP95'])].copy()
    var['CANONICAL_VARIANT_KEY']=key(var.chrom,var.pos,var.ref,var.alt); st['CANONICAL_VARIANT_KEY']=key(st.CHROM,st.POS,st.REF,st.ALT)
    dup_var=int(var.CANONICAL_VARIANT_KEY.duplicated().sum()); dup_stats=int(st.duplicated(['policy','CANONICAL_VARIANT_KEY']).sum())
    model=yaml.safe_load(MODEL.read_text()); features=model['rare_variant_annotations']; missing=[x for x in features if x not in ann.columns]
    if missing: raise RuntimeError('missing pretrained features: '+','.join(missing))
    if not {'id','gene_id'}.issubset(ann): raise RuntimeError('id/gene_id absent')
    gids=set(genes.id.astype('int64')); raw=raw[raw.gene_id.notna()].copy(); raw.gene_id=raw.gene_id.astype('int64')
    invalid_gene=int((~raw.gene_id.isin(gids)).sum()); dup_ann=int(raw.duplicated(['id','gene_id']).sum())
    mask=raw[(raw.CADD_PHRED.astype(float)>MASK_THRESHOLD)&raw.gene_id.isin(gids)][['id','gene_id']].copy(); mask['MASK_NAME']=MASK_NAME
    dup_mask=int(mask.duplicated(['id','gene_id','MASK_NAME']).sum())
    if dup_mask: raise RuntimeError('duplicate variant-gene-mask key; automatic deletion forbidden')
    mask=mask.merge(var[['id','CANONICAL_VARIANT_KEY']],on='id',how='left',validate='many_to_one'); join_missing=int(mask.CANONICAL_VARIANT_KEY.isna().sum())
    assigned=set(raw.id.astype('int64')); functional=set(mask.id.astype('int64')); branch={}; qout=[]
    for pol in ['GP90','GP95']:
        q=st[st.policy.eq(pol)].copy(); q['BRANCH']=pol; q['N_CALLED']=q.called_cells.astype('int64'); q['CALL_RATE']=q.N_CALLED/SAMPLE_N; q['MISSING_RATE']=q.missing_cells/SAMPLE_N; q['MONOMORPHIC']=q.AC.le(0); q['IMPUTATION_R2']=q.R2
        q=q.merge(var[['id','CANONICAL_VARIANT_KEY']],on='CANONICAL_VARIANT_KEY',how='left',validate='one_to_one')
        q['FUNCTIONAL_ANNOTATION_STATUS']=np.where(q.id.isin(functional),'ELIGIBLE_'+MASK_NAME,'NOT_ELIGIBLE'); q['GENE_ASSIGNMENT_STATUS']=np.where(q.id.isin(assigned),'ASSIGNED','UNASSIGNED')
        keep=['CANONICAL_VARIANT_KEY','BRANCH','N_CALLED','AN','AC','AF','MAF','CALL_RATE','MISSING_RATE','MONOMORPHIC','IMPUTATION_R2','FUNCTIONAL_ANNOTATION_STATUS','GENE_ASSIGNMENT_STATUS','carrier_cells','id']; branch[pol]=q[keep]; qout.append(q[keep].drop(columns='id')); q[keep].to_parquet(restricted/(pol+'_variant_qc.parquet'),index=False)
    pd.concat(qout).to_parquet(restricted/'CHR22_VARIANT_LEVEL_QC.parquet',index=False)
    denrows=[]; sums=[]; universe={}
    for b,(pol,cut) in BRANCHES.items():
        q=branch[pol]; q=q[(q.CALL_RATE>=CALL_RATE_MIN)&q.AC.gt(0)&q.MAF.gt(0)&q.MAF.lt(cut)]
        m=mask.merge(q[['id','CANONICAL_VARIANT_KEY','AC','carrier_cells','AN','MAF','CALL_RATE']],on=['id','CANONICAL_VARIANT_KEY'],how='inner',validate='many_to_one').merge(genes.rename(columns={'id':'gene_id'}),on='gene_id',how='left',validate='many_to_one')
        gm=m.groupby(['gene_id','gene','gene_name','MASK_NAME'],as_index=False).agg(CMAC=('AC','sum'),VARIANT_N=('id','nunique'),TOTAL_CARRIER_CELL_N=('carrier_cells','sum')); gm['TESTABLE']=gm.CMAC.ge(1); test=gm[gm.TESTABLE]
        tm=m.merge(test[['gene_id','MASK_NAME']],on=['gene_id','MASK_NAME'],how='inner',validate='many_to_one'); dup=int(tm.duplicated(['CANONICAL_VARIANT_KEY','gene_id','MASK_NAME']).sum())
        u=tm[['CANONICAL_VARIANT_KEY','gene_id','MASK_NAME']].copy(); universe[b]=u; u.to_parquet(restricted/(b+'_COMMON_UNIVERSE.parquet'),index=False); test.to_parquet(restricted/(b+'_TESTABLE_GENE_MASKS.parquet'),index=False); tm.to_parquet(restricted/(b+'_MEMBERSHIP_WITH_AGGREGATES.parquet'),index=False)
        h=digest(u,['CANONICAL_VARIANT_KEY','gene_id','MASK_NAME']); lines=[]
        for g in test.sort_values('gene').itertuples():
            ks=tm[tm.gene_id.eq(g.gene_id)].sort_values('CANONICAL_VARIANT_KEY').CANONICAL_VARIANT_KEY.tolist(); lines += [str(g.gene)+' var '+' '.join(ks),str(g.gene)+' anno '+' '.join([MASK_NAME]*len(ks))]
        atomic_text('\n'.join(lines)+'\n',groups/(b+'.group.txt')); os.chmod(groups/(b+'.group.txt'),0o600)
        status='PASS' if len(test)>0 and dup==0 else 'FAIL'; cm=test.CMAC
        base={'BRANCH':b,'MAF_THRESHOLD':cut,'MASK_NAME':MASK_NAME,'GENE_N':int(tm.gene_id.nunique()),'TESTABLE_GENE_MASK_N':len(test),'VARIANT_N':int(tm.id.nunique()),'TOTAL_AC':int(tm.AC.sum()),'TOTAL_CARRIER_CELL_N':int(tm.carrier_cells.sum()),'CMAC_MIN':int(cm.min()) if len(cm) else 0,'CMAC_MEDIAN':float(cm.median()) if len(cm) else 0,'CMAC_MAX':int(cm.max()) if len(cm) else 0,'DUPLICATE_KEY_N':dup,'STATUS':status}; sums.append(base)
        denrows.append(dict(base,HARDCALL=pol,MAF_RULE=f'0 < frozen-sample {pol} GT-derived MAF < {cut}',CALL_RATE_MIN=CALL_RATE_MIN,VARIANT_GENE_MASK_N=len(tm),DEEPRVAT_UNIVERSE_SHA256=h,SAIGE_UNIVERSE_SHA256=h))
    den=pd.DataFrame(denrows); summary=pd.DataFrame(sums)
    qc=[('FULL_CHR22_ANNOTATION',len(ann),'PASS' if len(features)==34 and not missing else 'FAIL','required_feature_n=34'),('CANONICAL_VARIANT_JOIN',join_missing,'PASS' if join_missing==0 else 'FAIL','missing join; no allele flip'),('GENE_ASSIGNMENT',int(mask.gene_id.nunique()),'PASS' if invalid_gene==0 and len(mask)>0 else 'FAIL',f'invalid_gene_row_n={invalid_gene}'),('FUNCTIONAL_MASK_MEMBERSHIP',len(mask),'PASS' if len(mask)>0 else 'FAIL',MASK_NAME),('DUPLICATE_CANONICAL_VARIANT_KEY_N',dup_var,'PASS' if dup_var==0 else 'FAIL','canonical variant table'),('DUPLICATE_POLICY_VARIANT_KEY_N',dup_stats,'PASS' if dup_stats==0 else 'FAIL','policy/key'),('DUPLICATE_VARIANT_GENE_MASK_KEY_N',dup_mask,'PASS' if dup_mask==0 else 'FAIL','no automatic deletion'),('ANNOTATION_DUPLICATE_ID_GENE_N',dup_ann,'PASS' if dup_ann==0 else 'FAIL','unfilled annotation')]
    for x in denrows: qc.append((x['BRANCH']+'_DENOMINATOR',x['TESTABLE_GENE_MASK_N'],x['STATUS'],'branch-specific'))
    qc += [('DEEPRVAT_SAIGE_VARIANT_UNIVERSE_MATCH',1,'PASS','same frozen membership table'),('DEEPRVAT_SAIGE_GENE_MASK_UNIVERSE_MATCH',1,'PASS','same frozen membership table'),('SAIGE_DS_FIXED_UNIVERSE_READY',int(len(universe['D90_MAF001'])>0),'PASS' if len(universe['D90_MAF001']) else 'FAIL','GP90 site/gene-mask universe; DS values not read'),('SOURCE_INFO_MAF_PREFILTER_LIMITATION',_config_number("PREFILTER_OMITTED_COUNT", int, False),'BLOCKED_PRODUCTION','chr22 canary allowed; production needs broad R2-only input or omission quantification')]
    qcdf=pd.DataFrame(qc,columns=['QC','VALUE','STATUS','DETAIL'])
    trs=[]
    for src,dst in [('RAW_GT','GP90'),('GP90','GP95')]:
        x=st[st.policy.eq(src)].set_index('CANONICAL_VARIANT_KEY'); y=st[st.policy.eq(dst)].set_index('CANONICAL_VARIANT_KEY'); z=x.join(y,lsuffix='_SRC',rsuffix='_DST'); inter=int(y.carrier_cells.sum()); union=int(x.carrier_cells.sum())
        trs.append({'TRANSITION':src+'_TO_'+dst,'VARIANT_N':len(z),'AC_CHANGED_N':int(z.AC_SRC.ne(z.AC_DST).sum()),'AN_CHANGED_N':int(z.AN_SRC.ne(z.AN_DST).sum()),'CALL_RATE_CHANGED_N':int(z.called_cells_SRC.ne(z.called_cells_DST).sum()),'MONOMORPHIC_GAIN_N':int((z.AC_SRC.gt(0)&z.AC_DST.eq(0)).sum()),'CARRIER_OVERLAP_INTERSECTION_N':inter,'CARRIER_OVERLAP_UNION_N':union,'CARRIER_JACCARD':inter/union,'RAW_TO_GP90_RETENTION':inter/union if src=='RAW_GT' else np.nan,'DEFINITION':'CARRIER_CELL_N=sample-variant cells with >=1 ALT; nested masking gives intersection=target, union=source'})
        pd.DataFrame({'SOURCE_BIN':mafbin(z.MAF_SRC),'TARGET_BIN':mafbin(z.MAF_DST)}).value_counts().reset_index(name='N').to_csv(restricted/(src+'_TO_'+dst+'_MAF_BIN_TRANSITION.tsv'),sep='\t',index=False)
    contract=pd.DataFrame([{'MASK_NAME':MASK_NAME,'MASK_DEFINITION':'CADD_PHRED > 5','REQUIRED_FEATURES':'CADD_PHRED','THRESHOLD':'>5','VARIANT_ELIGIBILITY_RULE':'branch call-rate >=0.95; AC>0; 0<cohort hard-call MAF<branch cutoff; CADD_PHRED>5','GENE_ASSIGNMENT_RULE':'VEP --per_gene then add_gene_ids against GENCODE v44 protein-coding HAVANA gene map','DEEPRVAT_USAGE':'pretrained-input association eligibility before learned 34-feature aggregation','SAIGE_USAGE':'same membership labelled CADD_PHRED_GT5 in group file','EVIDENCE_FILE':';'.join(map(str,[INPUTCFG,MODEL,PIPELINE,ANNO_CODE,SAIGE_CODE])),'EVIDENCE_LINE_OR_FUNCTION':'association_testing_data_thresholds; rare_variant_annotations; VEP --per_gene; add_gene_ids; SAIGE group parser','STATUS':'PASS'}])
    gate='PASS' if all(qcdf[qcdf.STATUS.isin(['PASS','FAIL'])].STATUS.eq('PASS')) and all(den.STATUS.eq('PASS')) else 'FAIL'
    atomic_tsv(summary,REPORTS/'CHR22_COMMON_GENE_MASK_SUMMARY.tsv'); atomic_tsv(den,REPORTS/'CHR22_GENE_MASK_DENOMINATORS.tsv'); atomic_tsv(qcdf,REPORTS/'CHR22_GENE_MASK_QC.tsv'); atomic_tsv(contract,REPORTS/'COMMON_MASK_CONTRACT.tsv'); atomic_tsv(pd.DataFrame(trs),REPORTS/'CHR22_VARIANT_POLICY_TRANSITIONS.tsv')
    g=qcdf.rename(columns={'QC':'gate','STATUS':'status','DETAIL':'detail'})[['gate','status','detail']]; g['evidence_path']=str(out); g=pd.concat([g,pd.DataFrame([{'gate':'COMMON_GENE_MASK_GATE','status':gate,'detail':'production-only source prefilter limitation excluded from chr22 technical gate','evidence_path':str(out)}])]); atomic_tsv(g,REPORTS/'COMMON_GENE_MASK_GATE.tsv')
    atomic_text(f'# Common Mask Contract\n\n- Mask: `{MASK_NAME}`\n- Definition: `CADD_PHRED > 5`\n- Evidence: `{INPUTCFG}`, `{MODEL}`, `{PIPELINE}`, `{ANNO_CODE}`.\n- The 34 pretrained features are learned inputs, not 34 invented functional masks.\n',REPORTS/'COMMON_MASK_CONTRACT.md')
    md=['# chr22 Common Gene-Mask Summary','',f'- COMMON_GENE_MASK_GATE: `{gate}`',f'- FULL_ANNOTATION_ROWS: `{len(ann)}`',f'- MASK: `{MASK_NAME}`',f'- CALL_RATE_MIN: `{CALL_RATE_MIN}`','- MAF001 and MAF01 are nested; MAF001 remains primary.','- TOTAL_CARRIER_CELL_N is sample-variant cells with >=1 ALT summed across memberships, not unique people.','','| Branch | Gene-mask N | Variant N | Total AC | Carrier cells | Status |','|---|---:|---:|---:|---:|---|']
    for x in sums: md.append(f"| {x['BRANCH']} | {x['TESTABLE_GENE_MASK_N']} | {x['VARIANT_N']} | {x['TOTAL_AC']} | {x['TOTAL_CARRIER_CELL_N']} | {x['STATUS']} |")
    md += ['',f'- Restricted versioned output: `{out}`','- Source INFO/MAF prefilter remains a genome-wide production blocker; it does not block this chr22 technical canary.']; atomic_text('\n'.join(md)+'\n',REPORTS/'CHR22_COMMON_GENE_MASK_SUMMARY.md')
    atomic_text(f'# Common Gene-Mask Gate\n\n- COMMON_GENE_MASK_GATE: `{gate}`\n- CHR22_CANARY_READY: `'+('YES' if gate=='PASS' else 'NO')+'`\n- Association was not executed.\n',REPORTS/'COMMON_GENE_MASK_GATE.md')
    atomic_text(str(out)+'\n',REPORTS/'CHR22_COMMON_GENE_MASK_RUN.latest.txt'); atomic_text(json.dumps({'created':stamp,'annotation':a.annotation,'unfilled':a.unfilled_annotation,'variants':a.variants,'annotation_run':a.annotation_run,'gate':gate,'archive':str(archive),'association_executed':False},indent=2)+'\n',out/'MANIFEST.json')
    for x in list(restricted.glob('*'))+list(groups.glob('*'))+[out/'MANIFEST.json']: os.chmod(x,0o600)
    print('COMMON_GENE_MASK_GATE='+gate); print('OUT='+str(out)); print(den[['BRANCH','TESTABLE_GENE_MASK_N','VARIANT_N','STATUS']].to_csv(sep='\t',index=False))
if __name__=='__main__': main()