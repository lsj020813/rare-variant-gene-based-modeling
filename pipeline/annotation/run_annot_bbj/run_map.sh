#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
cd ${PROJECT_ROOT}/work/run_annot_bbj
export TMPDIR=${PROJECT_ROOT}/work/tmp TMP=${PROJECT_ROOT}/work/tmp TEMP=${PROJECT_ROOT}/work/tmp
PY=python3
( ulimit -v 16000000; $PY map_bbj.py > map_full.log 2>&1; echo "map_full rc=$?" >> stage1.log ) &
sleep 3
bash ${PROJECT_ROOT}/work/run_band15/wd5.sh 'map_bb[j].py|pip_bbj_[A].py' ${PROJECT_ROOT}/work/run_annot_bbj 38 30 >> wd_map.log 2>&1 &
wait
