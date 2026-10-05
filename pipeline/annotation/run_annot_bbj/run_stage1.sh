#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
cd ${PROJECT_ROOT}/work/run_annot_bbj
export TMPDIR=${PROJECT_ROOT}/work/tmp TMP=${PROJECT_ROOT}/work/tmp TEMP=${PROJECT_ROOT}/work/tmp
PY=python3
bash ${PROJECT_ROOT}/work/run_band15/wd5.sh 'map_bb[j].py|pip_bbj_[A].py' ${PROJECT_ROOT}/work/run_annot_bbj 38 30 >> wd_stage1.log 2>&1 &
sleep 2
( ulimit -v 16000000; $PY map_bbj.py --smoke > map_smoke.log 2>&1; echo "map rc=$?" >> stage1.log ) &
( ulimit -v 16000000; $PY pip_bbj_A.py > pipA.log 2>&1; echo "pipA rc=$?" >> stage1.log ) &
wait
echo STAGE1_DONE >> stage1.log
