#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
cd ${PROJECT_ROOT}/work/run_annot_bbj
PY=python3
$PY pip_bbj_B.py --annot bbj_annot_sub.tsv.gz --outdir bbj_pip_sub --traits-out traits_sub.tsv > pipB_sub2.log 2>&1 || { echo "pipB_sub2 FAIL" >> stage5.log; exit 1; }
echo "pipB_sub2 ok" >> stage5.log
rm -rf ${PROJECT_ROOT}/work/run_band15/model_v10_out/e6_runs/smoke_sub_f0
bash ${PROJECT_ROOT}/work/run_annot_bbj/v10_run.sh smoke_sub_f0 _sub --fold 0 --phi-columns v9-eighteen --c1-ecdf-reference training --train-bands all
