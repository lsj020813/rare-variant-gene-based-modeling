#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
cd ${PROJECT_ROOT}/work/run_band15
KG=${PROJECT_ROOT}/work/ref15/annot/l1/known_genes_truthB.tsv
SP=$1
echo "[chain2] $(date '+%m-%d %H:%M') waiting spline pid $SP"
while kill -0 $SP 2>/dev/null; do sleep 60; done
W=$(pgrep -f "wd5[.]sh l1_trai" | head -1); [ -n "$W" ] && { kill -TERM $W 2>/dev/null; sleep 2; kill -9 $W 2>/dev/null; }
echo "[chain2] $(date '+%m-%d %H:%M') spline end: $(grep -c L1_DONE l1v8_spline_fold0.log) done, $(grep -c 'GATE FAIL\|Traceback' l1v8_spline_fold0.log) fail"
sleep 30
for arm in linear cluster nn; do
  echo "[chain2] $(date '+%m-%d %H:%M') start $arm"
  setsid nohup bash rt.sh l1v8_${arm}_fold0.log --fold 0 --arm $arm --threads 8 --known-genes $KG --inflation-range 0.8 1.2 --clamp-max 0.10 --corr-max 0.30 > /dev/null 2>&1 < /dev/null &
  sleep 90
  TP=$(pgrep -f "l1_trai[n]_v8.py --fold 0 --arm ${arm} " | head -1)
  if [ -z "$TP" ]; then echo "[chain2] $arm did not start (see log)"; sleep 10; continue; fi
  setsid nohup bash wd5.sh "l1_trai[n]_v8.py --fold" ${PROJECT_ROOT}/work/run_band15 38 170 >> wd_v8.log 2>&1 < /dev/null &
  sleep 5; WPID=$(pgrep -f "wd5[.]sh l1_trai" | head -1)
  echo "[chain2] $arm trainer=$TP watchdog=$WPID"
  while kill -0 $TP 2>/dev/null; do sleep 60; done
  [ -n "$WPID" ] && { kill -TERM $WPID 2>/dev/null; sleep 2; kill -9 $WPID 2>/dev/null; }
  echo "[chain2] $(date '+%m-%d %H:%M') end $arm: $(grep -c L1_DONE l1v8_${arm}_fold0.log) done, $(grep -c 'GATE FAIL\|Traceback' l1v8_${arm}_fold0.log) fail"
  sleep 30
done
echo "L1V8CHAIN_DONE"
