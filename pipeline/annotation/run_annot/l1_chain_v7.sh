#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
cd ${PROJECT_ROOT}/work/run_annot
export OMP_NUM_THREADS=8
PY=python3
for arm in spline nn linear; do
  echo "[chain] $(date '+%m-%d %H:%M') start $arm"
  $PY l1_train_v7.py --fold 0 --arm $arm --maxiter 60 --inner-maxiter 30 --null-perm 100 > l1v7_${arm}_fold0.log 2>&1 &
  TP=$!; sleep 5
  setsid nohup bash wd5.sh 'l1_train_v[7]' ${PROJECT_ROOT}/work/run_annot 38 130 >> wd5_launch.log 2>&1 </dev/null &
  WD=$!
  wait $TP; rc=$?
  kill -TERM $WD 2>/dev/null; sleep 1
  echo "[chain] $(date '+%m-%d %H:%M') end $arm rc=$rc $(grep -c L1_DONE l1v7_${arm}_fold0.log)"
done
echo L1V7CHAIN_DONE
