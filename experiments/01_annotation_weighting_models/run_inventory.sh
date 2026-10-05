#!/usr/bin/env bash
: "${CONDA_INIT_SCRIPT:?Set CONDA_INIT_SCRIPT}"
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -o pipefail
W=${PROJECT_ROOT}/work/run_l3
export TMPDIR=${PROJECT_ROOT}/work/tmp TMP=${PROJECT_ROOT}/work/tmp TEMP=${PROJECT_ROOT}/work/tmp
mkdir -p "$TMPDIR" "$W/out"
source ${CONDA_INIT_SCRIPT} && conda activate "${CONDA_ENV:?Set CONDA_ENV}"
cd "$W"
nice -n 19 ionice -c3 python3 l3_inventory.py \
  --root       ${PROJECT_ROOT}/work/ref \
  --alt-root   ${PROJECT_ROOT}/work/ref15 \
  --features-dir ${PROJECT_ROOT}/work/ref/features \
  --cadd-extract-dir ${PROJECT_ROOT}/work/ref/annot/extract \
  --resid-dir  ${PROJECT_ROOT}/work/ref/annot/resid \
  --truth-dir  ${PROJECT_ROOT}/work/ref/saige_step2_bwg \
  --gwas-dir   ${PROJECT_ROOT}/work/ref/gwas05 \
  --gencode    ${PROJECT_ROOT}/work/ref/deductive/gencode.sorted.gtf.gz \
  --band-vcf-dir ${PROJECT_ROOT}/work/ref/band_vcf \
  --phi-json   ${PROJECT_ROOT}/work/run_band15/model_v10_out/e6_runs/primary_f0/phi.json \
  --traits     tchl,htn,dm,lip \
  --threads    1 \
  --out        "$W/out"
rc=$?
echo "INVENTORY_RC=$rc"
exit $rc
