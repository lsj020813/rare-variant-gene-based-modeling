#!/usr/bin/env bash
: "${CONDA_INIT_SCRIPT:?Set CONDA_INIT_SCRIPT}"
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -o pipefail
W=${PROJECT_ROOT}/work/run_l3
export TMPDIR=${PROJECT_ROOT}/work/tmp TMP=${PROJECT_ROOT}/work/tmp TEMP=${PROJECT_ROOT}/work/tmp
mkdir -p "$TMPDIR" "$W/out"
source ${CONDA_INIT_SCRIPT} && conda activate "${CONDA_ENV:?Set CONDA_ENV}"
cd "$W"
COMMON="--root ${PROJECT_ROOT}/work/ref \
 --cadd-dir ${PROJECT_ROOT}/work/ref/features --cadd-source features \
 --chrom 22 --traits tchl --resid-dir ${PROJECT_ROOT}/work/ref/annot/resid \
 --ecdf-reference band-all --t1-partial-policy blank \
 --min-n 50000 --threads 4 --memory-gb 40 --min-avail-gb 80 \
 --resid-shuffle-seed 20260909 --perm-strata none"
nice -n 19 ionice -c3 python3 l3_burden.py $COMMON \
  --arms learned,flat,cadd_rank,pibar \
  --phi-json ${PROJECT_ROOT}/work/run_band15/model_v10_out/e6_runs/primary_f0/phi.json \
  --v10-script ${PROJECT_ROOT}/work/run_band15/model_v10_out/l1_train_v10.py \
  --pair-weight inverse-n-genes --perm-b 1 --perm-seed 777 --perm-block 1 \
  --pmode empirical --alpha 0.05 --multiple bonferroni-within-trait \
  --gene-batch 32 --out "$W/out" --tag smoke22_burden
rc1=$?; echo "SMOKE_BURDEN_RC=$rc1"
nice -n 19 ionice -c3 python3 l3_staaro.py $COMMON \
  --beta-weights '1,1;1,25' --acat-set annot-x-weight-x-test \
  --staar-null residual-analytic --skat-tail imhof \
  --perm-b 1 --perm-seed 777 --max-variants-per-gene 2000 \
  --out "$W/out" --tag smoke22_staaro
rc2=$?; echo "SMOKE_STAARO_RC=$rc2"
echo "SMOKE_DONE rc=$rc1/$rc2"
exit $(( rc1 | rc2 ))
