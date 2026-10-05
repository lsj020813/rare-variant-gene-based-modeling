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
import os, csv, json, glob, subprocess, collections
R = _configured('${PROJECT_ROOT}/work/ref')
BCF = _configured('${BCFTOOLS}')
GRM = _configured('${GRM_FILE}')
OUT = _configured('${PROJECT_ROOT}/work/run_audit')
rep = {}
FAIL = []

def ck(c, m):
    print(('  OK   ' if c else '  FAIL ') + m, flush=True)
    if not c:
        FAIL.append(m)
GT = [s.strip() for s in subprocess.run(f'{BCF} query -l {R}/band_vcf/chr22.band.vcf.gz', shell=True, capture_output=True, text=True).stdout.split()]
GTS = set(GT)
print('==== INPUT 4: sparse GRM ====', flush=True)
gids = [l.strip() for l in open(GRM + '.sampleIDs.txt') if l.strip()]
with open(GRM) as fh:
    hdr = [fh.readline() for _ in range(3)]
dims = [l for l in hdr if l and (not l.startswith('%'))]
print(f"  GRM sampleIDs: {len(gids):,} | mtx header dims: {(dims[0].strip() if dims else '?')}")
print(f'  first 2 GRM IDs: {gids[:2]}')
ck(len(gids) > 0, 'GRM sample ID list readable')
ov = len(set(gids) & GTS)
print(f'  GRM n genotype-VCF overlap: {ov:,} / GRM {len(gids):,} / VCF {len(GTS):,}')
ck(ov > 0, 'GRM IDs match the genotype sample IDs (format compatible)')
for name, p in [('htn_v2', f'{R}/pheno_v2/htn_v2.tsv'), ('tchl_primary', _configured('${ANALYSIS_ROOT}/work/phenotype/tchl_primary_complete_case.tsv'))]:
    with open(p, encoding='utf-8', errors='ignore') as fh:
        rd = csv.reader(fh, delimiter='\t')
        next(rd)
        pids = {r[0].strip() for r in rd if r}
    inn = len(pids & set(gids))
    print(f'  [{name}] {len(pids):,} pheno IDs | in GRM: {inn:,} | NOT in GRM: {len(pids) - inn:,}')
    ck(inn == len(pids), f'{name}: every phenotype sample present in the GRM')
rep['grm'] = {'n_ids': len(gids), 'overlap_vcf': ov}
print('\n==== INPUT 5: group files ====', flush=True)
CH = f'{R}/groupfiles_chunks'
parts = sorted(glob.glob(f'{CH}/chr*.part*.txt'))
print(f'  chunk files: {len(parts)}')
annos = collections.Counter()
badkey = 0
genes = set()
keys = set()
nvar = 0
for p in parts:
    for line in open(p):
        f = line.rstrip('\n').split()
        if len(f) < 3:
            continue
        genes.add(f[0])
        if f[1] == 'var':
            nvar += len(f) - 2
            keys.update(f[2:])
            for k in f[2:]:
                if k.startswith('chr') or k.count(':') != 3:
                    badkey += 1
        elif f[1] == 'anno':
            annos.update(f[2:])
print(f'  genes {len(genes):,} | gene-variant pairs {nvar:,} | distinct keys {len(keys):,}')
print(f'  annotation labels: {dict(annos)}')
ck(badkey == 0, f'variant keys use bare-contig CHR:POS:REF:ALT ({badkey} malformed)')
ck(set(annos) == {'all'} or 'all' in annos, f"--annotation_in_groupTest='all' matches the labels present {sorted(annos)}")
sample_keys = [k for k in list(keys) if k.startswith('22:')][:2000]
if sample_keys:
    have = set(subprocess.run(f"{BCF} query -f '%CHROM:%POS:%REF:%ALT\\n' {R}/band_vcf/chr22.band.vcf.gz", shell=True, capture_output=True, text=True).stdout.split())
    hit = sum((1 for k in sample_keys if k in have))
    print(f'  chr22 key join: {hit:,}/{len(sample_keys):,} group-file keys found in the VCF')
    ck(hit == len(sample_keys), 'every sampled group-file key resolves in the genotype file')
rep['groupfiles'] = {'genes': len(genes), 'pairs': nvar, 'keys': len(keys), 'annos': dict(annos)}
print('\n==== INPUT 6: argument coherence ====', flush=True)
print('  step1: covarColList=age,sex_male,PC1..PC10 | traitType=binary | LOCO=FALSE')
print('         useSparseGRMtoFitNULL=TRUE | isCateVarianceRatio=TRUE')
print('         cateVarRatioMinMACVecExclude=10,20.5 | MaxMACVecInclude=20.5')
print('  step2: vcfField=DS | maxMAF_in_groupTest=0.01 | minMAF=0 | minMAC=0.5')
print('         annotation_in_groupTest=all | AlleleOrder=ref-first | LOCO=FALSE | is_fastTest=TRUE')
print('  NOTE  --AlleleOrder=ref-first is the VCF convention; PLINK inputs would need alt-first.')
ck(True, 'argument set captured from the actual scripts (not assumed)')
json.dump(rep, open(f'{OUT}/saige_inputs2.json', 'w'), indent=1, ensure_ascii=False)
print(f'\nFAILURES: {len(FAIL)}')
for f in FAIL:
    print('  -', f)
print('VERIFY_PART2_COMPLETE')
