#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -u
cd ${PROJECT_ROOT}/work/run_band15
KG=${PROJECT_ROOT}/work/ref15/annot/l1/known_genes_truthB.tsv
PY=python3
D=${PROJECT_ROOT}/work/run_band15
F1=$1
TPID=""; WPID=""
cleanup() { [ -n "$WPID" ] && kill -TERM "$WPID" 2>/dev/null; [ -n "$TPID" ] && kill -TERM "$TPID" 2>/dev/null; sleep 2
            [ -n "$WPID" ] && kill -9 "$WPID" 2>/dev/null; [ -n "$TPID" ] && kill -9 "$TPID" 2>/dev/null
            echo "[chain4] $(date '+%m-%d %H:%M') cleanup done"; }
trap cleanup EXIT INT TERM

echo "[chain4] $(date '+%m-%d %H:%M') waiting fold 1 pid $F1 (not my child: exit code unavailable, judged by L1_DONE)"
while kill -0 $F1 2>/dev/null; do sleep 60; done
W1=$(pgrep -f "wd5[.]sh l1_trai\[n\]_v8.py --fold 1 --arm spline" | head -1); [ -n "$W1" ] && { kill -TERM $W1 2>/dev/null; sleep 2; kill -9 $W1 2>/dev/null; }
D1=$(grep -c '^L1_DONE' $D/l1v8_spline_fold1.log)
echo "[chain4] $(date '+%m-%d %H:%M') fold 1 ended: L1_DONE=$D1"
if [ "$D1" -eq 0 ]; then echo "[chain4] fold 1 technical failure — abort"; exit 13; fi
sleep 30

run_stage() {
  local LOG=$1 PAT=$2; shift 2
  local P=$D/tp_$LOG.pid; rm -f "$P"
  bash -c "export TMPDIR=${PROJECT_ROOT}/work/tmp; echo \$\$ > $P; cd $D/model_v8_out; exec $PY ./l1_train_v8.py $*" > $D/$LOG.log 2>&1 &
  TPID=$!
  for i in $(seq 1 30); do [ -s "$P" ] && break; sleep 1; done
  [ -s "$P" ] || { echo "[chain4] $LOG: no pid file — abort"; exit 11; }
  bash wd5.sh "$PAT" $D 38 170 >> $D/wd_v8.log 2>&1 &
  WPID=$!; sleep 8
  kill -0 "$WPID" 2>/dev/null || { echo "[chain4] $LOG: watchdog died — abort"; kill -TERM "$TPID"; sleep 3; kill -9 "$TPID" 2>/dev/null; exit 12; }
  echo "[chain4] $LOG trainer=$TPID watchdog=$WPID"
  wait "$TPID"; local RC=$?
  kill -TERM "$WPID" 2>/dev/null; sleep 2; kill -9 "$WPID" 2>/dev/null; WPID=""
  local DONE=$(grep -c '^L1_DONE' $D/$LOG.log)
  echo "[chain4] $(date '+%m-%d %H:%M') end $LOG: exit=$RC L1_DONE=$DONE"
  TPID=""
  if [ "$RC" -ne 0 ] || [ "$DONE" -eq 0 ]; then echo "[chain4] $LOG technical failure — abort"; exit 13; fi
  sleep 30
}

echo "[chain4] $(date '+%m-%d %H:%M') LC smoke (chr22, tiny) — gate before the real learning curve"
( cd $D/model_v8_out && export TMPDIR=${PROJECT_ROOT}/work/tmp && timeout 900 $PY ./l1_train_v8.py --smoke 22 --arm spline --learning-curve 0.5 --lam-divisors 10 --inner-maxiter 2 --inner-null-perm 5 --null-perm 5 --threads 4 ) > $D/lc_smoke.log 2>&1
if grep -q '^L1_DONE learning_curve' $D/lc_smoke.log; then
  echo "[chain4] LC smoke ok — start learning curve (fold 1 design)"
  run_stage l1v8_lc_fold1 "l1_trai[n]_v8.py --fold 1 --arm spline --learning-curve" \
  --fold 1 --arm spline --learning-curve 0.25 0.5 0.75 --lam-divisors 1 10 --inner-maxiter 60 --threads 8 --known-genes $KG
else
  echo "[chain4] LC smoke FAILED — skipping learning curve, continuing folds (see lc_smoke.log)"
fi

for k in 2 3 4; do
  echo "[chain4] $(date '+%m-%d %H:%M') start fold $k"
  run_stage l1v8_spline_fold$k "l1_trai[n]_v8.py --fold $k --arm spline" \
    --fold $k --arm spline --threads 8 --known-genes $KG --inflation-range 0.8 1.2 --clamp-max 0.10 --corr-max 0.30 \
    --lam-divisors 1 3.16227766 10 100 --inner-maxiter 60 --maxiter 150
done
echo "SPLINE_FOLDS_DONE"
