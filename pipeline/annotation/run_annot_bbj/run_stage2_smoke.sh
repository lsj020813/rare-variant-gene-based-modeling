#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
cd ${PROJECT_ROOT}/work/run_annot_bbj
export TMPDIR=${PROJECT_ROOT}/work/tmp TMP=${PROJECT_ROOT}/work/tmp TEMP=${PROJECT_ROOT}/work/tmp
PY=python3
( bash pool_bbj.sh "20" 1 > pool_smoke20.log 2>&1; echo "pool20 rc=$?" >> stage2.log ) &
( bash gt_pool_bbj.sh "20" 1 > gt_smoke20.log 2>&1; echo "gt20 rc=$?" >> stage2.log ) &
( ulimit -v 16000000; $PY rsq_bbj.py > rsq.log 2>&1; echo "rsq rc=$?" >> stage2.log ) &
sleep 3
bash ${PROJECT_ROOT}/work/run_band15/wd5.sh 'ext_bb[j].sh|gpn_bb[j].sh|t2_bb[j].py|rsq_bb[j].py|pip_bbj_[A].py' ${PROJECT_ROOT}/work/run_annot_bbj 38 30 >> wd_stage2.log 2>&1 &
wait
echo STAGE2_SMOKE_DONE >> stage2.log
