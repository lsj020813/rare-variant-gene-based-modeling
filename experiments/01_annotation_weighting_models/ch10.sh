#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -u
cd ${PROJECT_ROOT}/work/run_band15
KG=${PROJECT_ROOT}/work/ref15/annot/l1/known_genes_truthB.tsv
PY=python3
D=${PROJECT_ROOT}/work/run_band15
TPID=""; WPID=""
cleanup() {
  [ -n "$WPID" ] && kill -TERM "$WPID" 2>/dev/null
  [ -n "$TPID" ] && kill -TERM "$TPID" 2>/dev/null
  sleep 2
  [ -n "$WPID" ] && kill -9 "$WPID" 2>/dev/null
  [ -n "$TPID" ] && kill -9 "$TPID" 2>/dev/null
  echo "[chain3] $(date '+%m-%d %H:%M') cleanup done (own children only)"
}
trap cleanup EXIT INT TERM

for k in 1 2 3 4; do
  echo "[chain3] $(date '+%m-%d %H:%M') start fold $k"
  P=$D/tp_$k.pid; rm -f "$P"
  bash -c "export TMPDIR=${PROJECT_ROOT}/work/tmp; echo \$\$ > $P; cd $D/model_v8_out; \
    exec $PY ./l1_train_v8.py --fold $k --arm spline --threads 8 --known-genes $KG \
      --inflation-range 0.8 1.2 --clamp-max 0.10 --corr-max 0.30 \
      --lam-divisors 1 3.16227766 10 100 --inner-maxiter 60 --maxiter 150" \
    > $D/l1v8_spline_fold${k}.log 2>&1 &
  TPID=$!
  for i in $(seq 1 30); do [ -s "$P" ] && break; sleep 1; done
  RP=$(cat "$P" 2>/dev/null || echo "")
  if [ -z "$RP" ]; then echo "[chain3] fold $k: no pid file — abort chain"; exit 11; fi
  bash wd5.sh "l1_trai[n]_v8.py --fold $k --arm spline" $D 38 170 >> $D/wd_v8.log 2>&1 &
  WPID=$!
  sleep 8
  if ! kill -0 "$WPID" 2>/dev/null; then
    echo "[chain3] fold $k: watchdog failed to stay up — stopping trainer, abort chain"
    kill -TERM "$TPID" 2>/dev/null; sleep 3; kill -9 "$TPID" 2>/dev/null; exit 12
  fi
  echo "[chain3] fold $k trainer=$TPID(pidfile $RP) watchdog=$WPID"
  wait "$TPID"; RC=$?
  kill -TERM "$WPID" 2>/dev/null; sleep 2; kill -9 "$WPID" 2>/dev/null; WPID=""
  DONE=$(grep -c '^L1_DONE' $D/l1v8_spline_fold${k}.log)
  echo "[chain3] $(date '+%m-%d %H:%M') end fold $k: exit=$RC L1_DONE=$DONE"
  if [ "$RC" -ne 0 ] || [ "$DONE" -eq 0 ]; then
    echo "[chain3] fold $k technical failure (exit $RC) — abort chain per review REVIEW2 §4"
    TPID=""; exit 13
  fi
  TPID=""
  sleep 30
done
echo "SPLINE_FOLDS_DONE"
