#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
cd ${PROJECT_ROOT}/work/run_annot_bbj
export TMPDIR=${PROJECT_ROOT}/work/tmp TMP=${PROJECT_ROOT}/work/tmp TEMP=${PROJECT_ROOT}/work/tmp
PY=python3
( ulimit -v 16000000; $PY rsq_bbj.py > rsq.log 2>&1; echo "rsq rc=$?" >> stage3.log ) &
( ulimit -v 16000000; $PY pip_bbj_A.py > pipA2.log 2>&1; echo "pipA2 rc=$?" >> stage3.log ) &
sleep 3
bash ${PROJECT_ROOT}/work/run_band15/wd5.sh 'rsq_bb[j].py|pip_bbj_[A].py' ${PROJECT_ROOT}/work/run_annot_bbj 38 30 >> wd_stage3.log 2>&1 &
wait
echo STAGE3_DONE >> stage3.log
