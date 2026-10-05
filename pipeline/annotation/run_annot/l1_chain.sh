#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -u; cd ${PROJECT_ROOT}/work/run_annot
PY=python3
while [ "$(ls ${PROJECT_ROOT}/work/ref/annot/ds/*.done 2>/dev/null | wc -l)" -lt 22 ]; do sleep 120; done
echo "[chain] DS 22/22 at $(date '+%F %T')"; sleep 60
run(){ ARM=$1; F=$2
  [ -s ${PROJECT_ROOT}/work/ref/annot/l1/${ARM}_fold$F.json ] && { echo "[chain] ${ARM} fold$F skip"; return 0; }
  (setsid nohup bash wd5.sh 'l1_trai[n]' ${PROJECT_ROOT}/work/run_annot 38 160 >> wd5_launch.log 2>&1 </dev/null &)
  echo "[chain] ${ARM} fold$F start $(date '+%F %T')"
  export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
  /usr/bin/time -f "[${ARM} fold$F] rss=%MKB wall=%es" $PY l1_train.py --fold $F --arm $ARM --maxiter 60 > l1_${ARM}_fold$F.log 2>&1 || echo "[chain] ${ARM} fold$F FAIL"
  tail -2 l1_${ARM}_fold$F.log | cut -c1-200; sleep 20; }
for PASS in 1 2; do
  echo "[chain] pass $PASS"
  run nn 0; run spline 0; run linear 0
  for F in 1 2 3 4; do run nn $F; run spline $F; done
  for F in 1 2 3 4; do run linear $F; done
done
echo L1CHAIN_DONE
