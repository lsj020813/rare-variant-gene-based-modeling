#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
cd ${PROJECT_ROOT}/work/run_band15
KG=${PROJECT_ROOT}/work/ref15/annot/l1/known_genes_truthB.tsv
for k in 1 2 3 4; do
  echo "[chain3] $(date '+%m-%d %H:%M') start fold $k"
  setsid nohup bash rt.sh l1v8_spline_fold${k}.log --fold $k --arm spline --threads 8 --known-genes $KG \
    --inflation-range 0.8 1.2 --clamp-max 0.10 --corr-max 0.30 \
    --lam-divisors 1 3.16227766 10 100 --inner-maxiter 60 --maxiter 150 > /dev/null 2>&1 < /dev/null &
  sleep 120
  TP=$(pgrep -f "l1_trai[n]_v8.py --fold ${k} " | head -1)
  if [ -z "$TP" ]; then echo "[chain3] fold $k did not start"; sleep 10; continue; fi
  setsid nohup bash wd5.sh "l1_trai[n]_v8.py --fold" ${PROJECT_ROOT}/work/run_band15 38 170 >> wd_v8.log 2>&1 < /dev/null &
  sleep 5; WPID=$(pgrep -f "wd5[.]sh l1_trai" | head -1)
  echo "[chain3] fold $k trainer=$TP watchdog=$WPID"
  while kill -0 $TP 2>/dev/null; do sleep 120; done
  [ -n "$WPID" ] && { kill -TERM $WPID 2>/dev/null; sleep 2; kill -9 $WPID 2>/dev/null; }
  echo "[chain3] $(date '+%m-%d %H:%M') end fold $k: $(grep -c L1_DONE l1v8_spline_fold${k}.log) done, $(grep -c 'GATE FAIL\|Traceback' l1v8_spline_fold${k}.log) fail"
  sleep 30
done
echo "SPLINE_FOLDS_DONE"
