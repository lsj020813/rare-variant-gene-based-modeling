#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -u
cd ${PROJECT_ROOT}/work/run_band15
KG=${PROJECT_ROOT}/work/ref15/annot/l1/known_genes_truthB.tsv
PY=python3
D=${PROJECT_ROOT}/work/run_band15
TPID=""; WPID=""
cleanup() { [ -n "$WPID" ] && kill -TERM "$WPID" 2>/dev/null; [ -n "$TPID" ] && kill -TERM "$TPID" 2>/dev/null; sleep 2
            [ -n "$WPID" ] && kill -9 "$WPID" 2>/dev/null; [ -n "$TPID" ] && kill -9 "$TPID" 2>/dev/null
            echo "[chain5] $(date '+%m-%d %H:%M') cleanup done"; }
trap cleanup EXIT INT TERM

echo "[chain5] $(date '+%m-%d %H:%M') waiting for common05 22/22 done + no extractors"
for i in $(seq 1 1440); do
  n=$(ls ${PROJECT_ROOT}/work/ref/common05/*.done 2>/dev/null | wc -l)
  b=$(pgrep -c -f 'bcftools view -i INFO/MAF' || true)
  [ "$n" -ge 22 ] && [ "${b:-0}" -eq 0 ] && break
  sleep 60
done
[ "$n" -ge 22 ] || { echo "[chain5] ABORT: extraction not finished after 24h (done=$n)"; exit 2; }
echo "[chain5] $(date '+%m-%d %H:%M') extraction complete (done=$n) — start folds"
sleep 30

run_stage() {
  local LOG=$1 PAT=$2; shift 2
  local P=$D/tp_$LOG.pid; rm -f "$P"
  bash -c "export TMPDIR=${PROJECT_ROOT}/work/tmp; echo \$\$ > $P; cd $D/model_v8_out; exec $PY ./l1_train_v8.py $*" > $D/$LOG.log 2>&1 &
  TPID=$!
  for i in $(seq 1 30); do [ -s "$P" ] && break; sleep 1; done
  [ -s "$P" ] || { echo "[chain5] $LOG: no pid file — abort"; exit 11; }
  bash wd5.sh "$PAT" $D 38 170 >> $D/wd_v8.log 2>&1 &
  WPID=$!; sleep 8
  kill -0 "$WPID" 2>/dev/null || { echo "[chain5] $LOG: watchdog died — abort"; kill -TERM "$TPID"; sleep 3; kill -9 "$TPID" 2>/dev/null; exit 12; }
  echo "[chain5] $LOG trainer=$TPID watchdog=$WPID"
  wait "$TPID"; local RC=$?
  kill -TERM "$WPID" 2>/dev/null; sleep 2; kill -9 "$WPID" 2>/dev/null; WPID=""
  local DONE=$(grep -c '^L1_DONE' $D/$LOG.log)
  echo "[chain5] $(date '+%m-%d %H:%M') end $LOG: exit=$RC L1_DONE=$DONE"
  TPID=""
  if [ "$RC" -ne 0 ] || [ "$DONE" -eq 0 ]; then echo "[chain5] $LOG technical failure — abort"; exit 13; fi
  sleep 30
}

for k in 3 4; do
  echo "[chain5] $(date '+%m-%d %H:%M') start fold $k"
  run_stage l1v8_spline_fold$k "l1_trai[n]_v8.py --fold $k --arm spline" \
    --fold $k --arm spline --threads 8 --known-genes $KG --inflation-range 0.8 1.2 --clamp-max 0.10 --corr-max 0.30 \
    --lam-divisors 1 3.16227766 10 100 --inner-maxiter 60 --maxiter 150
done
echo "SPLINE_FOLDS_DONE"
