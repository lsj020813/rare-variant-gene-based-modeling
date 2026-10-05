#!/usr/bin/env bash
: "${CONDA_INIT_SCRIPT:?Set CONDA_INIT_SCRIPT}"
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -o pipefail
W=${PROJECT_ROOT}/work/run_l3
export TMPDIR=${PROJECT_ROOT}/work/tmp TMP=$TMPDIR TEMP=$TMPDIR
mkdir -p "$TMPDIR" "$W/out"
source ${CONDA_INIT_SCRIPT} && conda activate "${CONDA_ENV:?Set CONDA_ENV}"
cd "$W"
nice -n 19 ionice -c3 python3 l3_burden.py \
  --root ${PROJECT_ROOT}/work/ref \
  --cadd-dir ${PROJECT_ROOT}/work/ref/features --cadd-source features \
  --chrom 19 --traits tchl --resid-dir ${PROJECT_ROOT}/work/ref/annot/resid \
  --arms learned,flat,cadd_rank,pibar \
  --phi-json ${PROJECT_ROOT}/work/run_band15/model_v10_out/e6_runs/primary_f0/phi.json \
  --v10-script ${PROJECT_ROOT}/work/run_band15/model_v10_out/l1_train_v10.py \
  --pair-weight inverse-n-genes \
  --ecdf-reference band-all --t1-partial-policy blank \
  --perm-b 20 --perm-seed 20260909 --perm-block 10 \
  --pmode moment --alpha 0.05 --multiple none \
  --min-n 50000 --gene-batch 32 --threads 4 --memory-gb 40 --min-avail-gb 80 \
  --resid-shuffle-seed none --perm-strata none --allow-real-residuals \
  --out "$W/out" --tag preview19_burden
rc=$?; echo "PREVIEW_BURDEN_RC=$rc"
[ $rc -ne 0 ] && exit $rc
nice -n 19 ionice -c3 python3 l3_metrics.py \
  --gene-trait out/preview19_burden.gene_trait.tsv \
  --burden-summary out/preview19_burden.SUMMARY.json \
  --staaro none \
  --truth-dir ${PROJECT_ROOT}/work/ref/saige_step2_bwg \
  --truth-group all --truth-maxmaf 0.01 --truth-alpha 2.5e-6 \
  --reference-arm flat --p-column p \
  --alpha-grid 0.05,0.01,0.001,0.0001,2.5e-6 --topk-grid 10,25,50,100 \
  --out "$W/out" --tag preview19_metrics
echo "PREVIEW_METRICS_RC=$?"
