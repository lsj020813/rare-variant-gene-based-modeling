#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -euo pipefail
BCF=${BCFTOOLS:-bcftools}
BAND=${PROJECT_ROOT}/work/ref/band_vcf/chr21.band.vcf.gz
GF=${PROJECT_ROOT}/work/ref/groupfiles_v1/chr21.groupfile_V1.txt
PILOT=${ANALYSIS_ROOT:?Set ANALYSIS_ROOT}
MODEL=$PILOT/work/saige_gene/grch38_chr22/step1/tchl_primary.rda
VR=$PILOT/work/saige_gene/grch38_chr22/step1/tchl_primary_with_vr_chr22_gp90_retry3.varianceRatio.txt
OUT=${PROJECT_ROOT}/work/ref/b9_smoke
mkdir -p "$OUT"
export UDOCKER_DIR=${CONTAINER_STATE_DIR:?Set CONTAINER_STATE_DIR}

python3 - <<'PYEOF'
import os as _config_os
import re as _config_re

def _required(name):
    value = _config_os.environ.get(name)
    if value is None or not value.strip():
        raise RuntimeError(f"Set {name} before running this script")
    return value

def _required_int(name):
    value = int(_required(name))
    if value <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return value

def _configured(value):
    defaults = {'BCFTOOLS': 'bcftools', 'PLINK2': 'plink2', 'PYTHON': 'python3', 'TABIX': 'tabix', 'BGZIP': 'bgzip', 'SAMTOOLS': 'samtools', 'BEDTOOLS': 'bedtools', 'CROSSMAP': 'CrossMap', 'UDOCKER_BIN': 'udocker'}
    return _config_re.sub(r"\$\{([A-Z][A-Z0-9_]*)\}|\$([A-Z][A-Z0-9_]*)", lambda m: _config_os.environ.get(m.group(1) or m.group(2)) or defaults.get(m.group(1) or m.group(2)) or _required(m.group(1) or m.group(2)), value)

def _source_path(base, value):
    return value if _config_os.path.isabs(value) else _config_os.path.join(base, value)


import gzip, math, subprocess
GF = _configured('${PROJECT_ROOT}/work/ref/groupfiles_v1/chr21.groupfile_V1.txt')
OUT = _configured('${PROJECT_ROOT}/work/ref/b9_smoke')
lines = open(GF).read().splitlines()
genes = {}
for i in range(0, len(lines), 2):
    g = lines[i].split(' ', 1)[0]
    genes[g] = (lines[i], lines[i + 1])
sel = sorted(genes)[:30]
# MAF lookup for all variants in selected genes
keys = set()
for g in sel:
    keys.update(genes[g][0].split()[2:])
import re
maf = {}
p = subprocess.Popen([_configured('${BCFTOOLS}'), 'query', '-f', '%CHROM:%POS:%REF:%ALT\t%INFO/MAF\n', _configured('${PROJECT_ROOT}/work/ref/band_vcf/chr21.band.vcf.gz')], stdout=subprocess.PIPE, text=True)
for line in p.stdout:
    k, m = line.rstrip('\n').split('\t')
    if k in keys:
        maf[k] = float(m)
p.stdout.close()
p.wait()
missing = [k for g in sel for k in genes[g][0].split()[2:] if k not in maf]
assert not missing, f'{len(missing)} keys missing MAF'
from math import lgamma

def beta_w(m, a=1.0, b=25.0):
    # dbeta(maf, 1, 25) as SAIGE uses
    return math.exp(lgamma(a + b) - lgamma(a) - lgamma(b) + (a - 1) * math.log(max(m, 1e-12)) + (b - 1) * math.log(max(1 - m, 1e-12)))
with open(f'{OUT}/gf_R1.txt', 'w') as f1, open(f'{OUT}/gf_R2.txt', 'w') as f2, open(f'{OUT}/gf_R3.txt', 'w') as f3:
    for g in sel:
        var, anno = genes[g]
        ks = var.split()[2:]
        for f in (f1, f2, f3):
            f.write(var + '\n')
            f.write(anno + '\n')
        f2.write(f'{g} weight ' + ' '.join(['1.0'] * len(ks)) + '\n')
        f3.write(f'{g} weight ' + ' '.join((f'{beta_w(maf[k]):.6g}' for k in ks)) + '\n')
print('genes:', len(sel), '| variants:', sum((len(genes[g][0].split()) - 2 for g in sel)))
PYEOF

run(){
  local TAG=$1 GFX=$2 EXTRA=${3:-}
  local O=$OUT/smoke_$TAG
  PROOT_NO_SECCOMP=1 nice -n 15 ${UDOCKER_BIN:-udocker} run --volume=/data:/data ${SAIGE_IMAGE:?Set SAIGE_IMAGE} \
    step2_SPAtests.R \
    --vcfFile=$BAND --vcfFileIndex=$BAND.csi --vcfField=DS --chrom=21 \
    --minMAF=0 --minMAC=0.5 --maxMAF_in_groupTest=0.01 \
    --GMMATmodelFile=$MODEL --varianceRatioFile=$VR \
    --SAIGEOutputFile=$O --groupFile=$GFX \
    --annotation_in_groupTest="all" $EXTRA \
    --is_fastTest=TRUE --LOCO=FALSE --AlleleOrder=ref-first > $O.log 2>&1 \
    || { echo "FAIL $TAG"; tail -3 $O.log; return 5; }
  local NL=$(wc -l < $O); [ "$NL" -ge 3 ] || { echo "FAIL $TAG: $NL rows"; return 6; }
  echo "[${TAG}] OK rows=$NL"
}
run R1 $GF &
run R2 $OUT/gf_R2.txt &
run R3 $OUT/gf_R3.txt &
run R4 $GF "--is_no_weight_in_groupTest=TRUE" &
wait
echo B9_SMOKE_COMPLETE
