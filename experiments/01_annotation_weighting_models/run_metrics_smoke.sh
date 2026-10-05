#!/usr/bin/env bash
: "${CONDA_INIT_SCRIPT:?Set CONDA_INIT_SCRIPT}"
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -o pipefail
W=${PROJECT_ROOT}/work/run_l3
export TMPDIR=${PROJECT_ROOT}/work/tmp TMP=$TMPDIR TEMP=$TMPDIR
source ${CONDA_INIT_SCRIPT} && conda activate "${CONDA_ENV:?Set CONDA_ENV}"
cd "$W"
nice -n 19 ionice -c3 python3 l3_metrics.py \
  --gene-trait out/smoke22_burden.gene_trait.tsv \
  --burden-summary out/smoke22_burden.SUMMARY.json \
  --staaro out/smoke22_staaro.staaro.tsv \
  --truth-dir ${PROJECT_ROOT}/work/ref/saige_step2_bwg \
  --truth-group all --truth-maxmaf 0.01 --truth-alpha 2.5e-6 \
  --reference-arm flat --p-column p \
  --alpha-grid 0.05,0.01,0.001,0.0001,2.5e-6 --topk-grid 10,25,50,100 \
  --out "$W/out" --tag smoke22_metrics
echo "METRICS_RC=$?"
