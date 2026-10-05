#!/usr/bin/env python3
import os as _config_os
import re as _config_re

def _required(name):
    value = _config_os.environ.get(name)
    if value is None or not value.strip():
        raise RuntimeError(f'Set {name} before running this script')
    return value

def _required_int(name):
    value = int(_required(name))
    if value <= 0:
        raise ValueError(f'{name} must be a positive integer')
    return value

def _configured(value):
    defaults = {'BCFTOOLS': 'bcftools', 'PLINK2': 'plink2', 'PYTHON': 'python3', 'TABIX': 'tabix', 'BGZIP': 'bgzip', 'SAMTOOLS': 'samtools', 'BEDTOOLS': 'bedtools', 'CROSSMAP': 'CrossMap', 'UDOCKER_BIN': 'udocker'}
    def resolve(match):
        name = match.group(1) or match.group(2)
        if name in _config_os.environ:
            return _required(name)
        return defaults.get(name) or _required(name)
    return _config_re.sub(r"\$\{([A-Z][A-Z0-9_]*)\}|\$([A-Z][A-Z0-9_]*)", resolve, value)

def _source_path(base, value):
    return value if _config_os.path.isabs(value) else _config_os.path.join(base, value)
import os, sys, csv, json, glob, gzip, subprocess, collections, statistics as st
R = _configured('${PROJECT_ROOT}/work/ref')
B = _configured('${PHENO_DIR}')
BCF = _configured('${BCFTOOLS}')
OUT = _configured('${PROJECT_ROOT}/work/run_audit')
os.makedirs(OUT, exist_ok=True)
rep, FAIL = ({}, [])

def hd(msg):
    print(f"\n{'=' * 4} {msg} {'=' * 4}", flush=True)

def ck(cond, msg):
    print(('  OK   ' if cond else '  FAIL ') + msg, flush=True)
    if not cond:
        FAIL.append(msg)
    return cond
hd('INPUT 1: genotype (--vcfFile, --vcfFileIndex, --vcfField, --chrom)')
g = {}
samp = None
for n in range(1, 23):
    v = f'{R}/band_vcf/chr{n}.band.vcf.gz'
    i = v + '.csi'
    cnt = subprocess.run(f'{BCF} index -n {v}', shell=True, capture_output=True, text=True).stdout.strip()
    g[f'chr{n}'] = {'vcf': os.path.exists(v), 'csi': os.path.exists(i), 'n_records': int(cnt) if cnt.isdigit() else None}
missing = [k for k, v in g.items() if not (v['vcf'] and v['csi'])]
ck(not missing, f'all 22 band VCFs + .csi present (missing: {missing})')
tot = sum((v['n_records'] or 0 for v in g.values()))
print(f'  total band records: {tot:,}')
s1 = subprocess.run(f'{BCF} query -l {R}/band_vcf/chr1.band.vcf.gz', shell=True, capture_output=True, text=True).stdout.split()
s22 = subprocess.run(f'{BCF} query -l {R}/band_vcf/chr22.band.vcf.gz', shell=True, capture_output=True, text=True).stdout.split()
ck(s1 == s22, f'sample order identical chr1 vs chr22 ({len(s1):,} samples)')
GT = {x.split('_')[0] for x in s1}
rep['input1_genotype'] = {'per_chrom': g, 'total_records': tot, 'n_samples': len(s1)}
fmt = subprocess.run(f"{BCF} view -h {R}/band_vcf/chr22.band.vcf.gz | grep -c '^##FORMAT=<ID=DS'", shell=True, capture_output=True, text=True).stdout.strip()
ck(fmt == '1', f'FORMAT/DS declared in header (--vcfField=DS)')
row = subprocess.run(f"{BCF} query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%DS]\\n' {R}/band_vcf/chr22.band.vcf.gz 2>/dev/null | head -1", shell=True, capture_output=True, text=True).stdout.split('\t')
ck(len(row) > 4 and row[4].strip() not in ('', '.'), f"DS values readable (first: {(row[4].strip() if len(row) > 4 else 'NONE')})")
ck(not row[0].startswith('chr'), f"--chrom convention: contig is '{row[0]}' (bare, no 'chr')")
rep['input1_contig_style'] = row[0]
afs = subprocess.run(f"{BCF} query -f '%INFO/MAF\\n' {R}/band_vcf/chr22.band.vcf.gz 2>/dev/null | head -20000", shell=True, capture_output=True, text=True).stdout.split()
vals = [float(x) for x in afs if x not in ('', '.')]
if vals:
    print(f'  MAF in file: min {min(vals):.5f} max {max(vals):.5f} (n={len(vals):,} sampled)')
    ck(min(vals) >= 0.0009 and max(vals) <= 0.0101, f'MAF within declared band 0.001-0.01')
    rep['input1_maf'] = {'min': min(vals), 'max': max(vals)}
else:
    print('  MAF INFO tag absent — maxMAF_in_groupTest will use computed AF')
hd('INPUT 2: phenotype (--phenoFile, --phenoCol, --covarColList, --sampleIDColinphenoFile)')
ph = {}
for name, p in [('htn_v2', f'{R}/pheno_v2/htn_v2.tsv'), ('dm_v2', f'{R}/pheno_v2/dm_v2.tsv'), ('lip_v2', f'{R}/pheno_v2/lip_v2.tsv'), ('tchl_primary', _configured('${ANALYSIS_ROOT}/work/phenotype/tchl_primary_complete_case.tsv'))]:
    if not os.path.exists(p):
        ph[name] = {'exists': False}
        ck(False, f'{name}: file missing')
        continue
    with open(p, encoding='utf-8', errors='ignore') as fh:
        rd = csv.reader(fh, delimiter='\t')
        hdr = next(rd)
        rows = [r for r in rd if r]
    ids = {r[0].split('_')[0] for r in rows}
    ycol = 'y' if 'y' in hdr else 'TCHL_rint' if 'TCHL_rint' in hdr else None
    yi = hdr.index(ycol) if ycol else None
    yv = [r[yi] for r in rows if yi is not None and len(r) > yi]
    nmiss = sum((1 for v in yv if v.strip() in ('', 'NA', '.')))
    ph[name] = {'exists': True, 'n': len(rows), 'cols': hdr, 'ycol': ycol, 'n_in_genotype': len(ids & GT), 'n_missing_y': nmiss, 'id_format': rows[0][0]}
    print(f"  [{name}] n={len(rows):,} in-genotype={len(ids & GT):,} ycol={ycol} missing_y={nmiss} id='{rows[0][0]}'")
    print(f"      covars available: {[c for c in hdr if c not in ('sample_id', ycol)][:16]}")
    ck(len(ids & GT) == len(ids), f'{name}: every phenotype ID is present in the genotypes')
rep['input2_pheno'] = ph
sets = {n: [c for c in v['cols'] if c.startswith(('age', 'sex', 'PC'))] for n, v in ph.items() if v.get('exists')}
uniq = {tuple(v) for v in sets.values()}
ck(len(uniq) == 1, f'covariate columns identical across traits ({ {n: len(v) for n, v in sets.items()}})')
for n, v in sets.items():
    print(f'      {n}: {v}')
rep['input2_covars'] = sets
hd('INPUT 3: null model + variance ratio')
nm = {}
for t, rda, vr in [('htn', f'{R}/saige_step1/htn_v2.rda', f'{R}/saige_step1/htn_v2.varianceRatio.txt'), ('dm', f'{R}/saige_step1/dm_v2.rda', f'{R}/saige_step1/dm_v2.varianceRatio.txt'), ('lip', f'{R}/saige_step1/lip_v2.rda', f'{R}/saige_step1/lip_v2.varianceRatio.txt'), ('tchl', _configured('${ANALYSIS_ROOT}/work/saige_gene/grch38_chr22/step1/tchl_primary.rda'), _configured('${ANALYSIS_ROOT}/work/saige_gene/grch38_chr22/step1/tchl_primary_with_vr_chr22_gp90_retry3.varianceRatio.txt'))]:
    e1, e2 = (os.path.exists(rda), os.path.exists(vr))
    body = open(vr).read().strip() if e2 else ''
    nm[t] = {'rda': rda, 'rda_exists': e1, 'rda_mb': round(os.path.getsize(rda) / 1000000.0, 1) if e1 else None, 'vr': vr, 'vr_exists': e2, 'vr_lines': body.splitlines() if body else []}
    print(f"  [{t}] rda={('Y' if e1 else 'N')} ({nm[t]['rda_mb']}MB)  vr={('Y' if e2 else 'N')}")
    for l in nm[t]['vr_lines'][:5]:
        print(f'      {l}')
    ck(e1 and e2 and body, f'{t}: null model and variance ratio both present and non-empty')
rep['input3_nullmodel'] = nm
hd('INPUT 4: sparse GRM (--sparseGRMFile, --sparseGRMSampleIDFile)')
cands = glob.glob('/data/**/**sparseGRM*.mtx', recursive=False) + glob.glob(_configured('${DATA_ROOT}/GWAS/Saige/sparseGRM/*sparseGRM*'), recursive=False) + glob.glob(f'{R}/**/*sparseGRM*', recursive=True)
cands = sorted(set(cands))
for p in cands[:8]:
    print(f'  found: {p}  ({os.path.getsize(p) / 1000000.0:.1f}MB)')
rep['input4_grm_candidates'] = cands
ck(bool(cands), 'sparse GRM file located')
json.dump(rep, open(f'{OUT}/saige_inputs.json', 'w'), indent=1, ensure_ascii=False)
print(f"\n{'=' * 40}\nFAILURES: {len(FAIL)}")
for f in FAIL:
    print('  -', f)
print('VERIFY_PART1_COMPLETE')
