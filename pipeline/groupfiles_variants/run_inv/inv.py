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
import os, glob, json, subprocess
ROOT = _configured('${PROJECT_ROOT}')
ANN = _configured('${ANALYSIS_ROOT}/resources/annotation_data')
W = f'{ROOT}/work'
PLAN = {'A_position_element': [('ENCODE_cCRE', [f'{W}/tmp/ccre/GRCh38-cCREs.V4.bed']), ('Gnocchi', [f'{W}/tmp/gnocchi/gnocchi_1kb_qc.txt.gz']), ('ReMap_TF_summary', [f'{W}/ref/remap/*.bed*', f'{W}/tmp/remap/*']), ('CpG_island', [f'{W}/ref/cpg/*', f'{W}/tmp/cpg*/*'])], 'B_tissue_context': [('rE2G', [f'{W}/ref/re2g/*.bed.gz']), ('CATlas', [f'{W}/ref/catlas/*', f'{W}/tmp/catlas/*'])], 'C_effect_size': [('CADD_wg_SNV', [f'{ANN}/cadd/whole_genome_SNVs.tsv.gz']), ('CADD_indel', [f'{ANN}/cadd/gnomad.genomes.r3.0.indel.tsv.gz']), ('annotations_45col', [_configured('${ANALYSIS_ROOT}/**/annotations.parquet')]), ('AlphaMissense', [f'{ANN}/AlphaMissense/*.tsv.gz']), ('primateAI', [f'{ANN}/primateAI/*.bgz', f'{ANN}/primateAI/*.gz'])], 'D_effect_direction': [('SpliceAI_raw_snv', [f'{ANN}/spliceAI/spliceai_scores.raw.snv.hg38.vcf.gz']), ('SpliceAI_raw_indel', [f'{ANN}/spliceAI/spliceai_scores.raw.indel.hg38.vcf.gz']), ('AlphaGenome', [f'{W}/ref/alphagenome/*', f'{W}/tmp/alphagenome/*']), ('TF_motif_JASPAR', [f'{W}/ref/jaspar/*', f'{W}/tmp/jaspar/*', f'{W}/ref/motif/*']), ('pLoF_label_VEP', [f'{ANN}/gencode/*.gtf.gz']), ('uAUG_UTRannotator', [f'{W}/ref/utrannotator/*', f'{W}/tmp/utr*/*']), ('Enformer_SAD', [f'{W}/ref/enformer*', f'{W}/tmp/enformer*'])], 'E_conservation': [('GPN_MSA', [f'{W}/ref/gpnmsa/*', f'{W}/tmp/gpn*/*']), ('aPC_Conservation', [f'{W}/ref/apc/*', f'{W}/tmp/apc*/*']), ('phyloP', [f'{W}/ref/phylop/*', f'{W}/tmp/phylop/*']), ('phastCons', [f'{W}/ref/phastcons/*', f'{W}/tmp/phastcons/*'])], 'F_splicing': [('Pangolin', [f'{W}/ref/pangolin/*', f'{W}/tmp/pangolin/*']), ('AbSplice', [f'{ANN}/absplice/*'])], 'K_rna_level': [('DeepRipe', [f'{W}/ref/deepripe*', f'{W}/tmp/deepripe*', f'{ROOT}/.kipoi/models/DeepRipe*']), ('TargetScan', [f'{W}/ref/targetscan/*', f'{W}/tmp/targetscan/*'])], 'L_technical_gate': [('mappability', [f'{W}/ref/mappability/*', f'{W}/tmp/mappab*/*']), ('repeatmasker', [f'{W}/ref/repeats/*', f'{W}/tmp/repeat*/*']), ('mut_rate', [f'{W}/ref/mutrate/*', f'{W}/tmp/mutrate/*'])]}
out = {}
for topic, cards in PLAN.items():
    out[topic] = {}
    for name, pats in cards:
        found = []
        for p in pats:
            for f in glob.glob(p, recursive=True)[:40]:
                try:
                    sz = os.path.getsize(f)
                except OSError:
                    continue
                if os.path.isfile(f) and sz > 1000000:
                    found.append((f, sz))
        found.sort(key=lambda x: -x[1])
        out[topic][name] = {'status': 'HAVE' if found else 'MISSING', 'n_files': len(found), 'largest': found[0][0] if found else None, 'bytes': found[0][1] if found else 0}
json.dump(out, open(_configured('${PROJECT_ROOT}/work/ref/b6_inventory.json'), 'w'), indent=1)
have = miss = 0
for topic, cards in out.items():
    print(f'\n=== {topic} ===')
    for n, v in cards.items():
        mark = 'HAVE ' if v['status'] == 'HAVE' else 'MISS '
        sz = f"{v['bytes'] / 2 ** 30:.2f}G" if v['bytes'] else '-'
        print(f"  {mark}{n:<22}{sz:>9}  {(v['largest'] or '')[-58:]}")
        have += v['status'] == 'HAVE'
        miss += v['status'] == 'MISSING'
print(f'\nHAVE {have} | MISSING {miss}')
