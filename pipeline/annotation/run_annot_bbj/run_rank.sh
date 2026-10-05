#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
cd ${PROJECT_ROOT}/work/run_annot_bbj
PY=python3
export TMPDIR=${PROJECT_ROOT}/work/run_band15/model_v10_out/.scratch_e6/rank_primary_f0; mkdir -p $TMPDIR
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
( $PY rank_metrics.py primary_f0 > rank_primary_f0.log 2>&1; echo "rank rc=$?" >> rank_primary_f0.log ) &
sleep 5
bash ${PROJECT_ROOT}/work/run_annot_bbj/wd_rss.sh 'rank_metric[s].py primary_f0' ${PROJECT_ROOT}/work/run_band15/model_v10_out/e6_runs/primary_f0/wd_rank.log 60 6 &
wait
