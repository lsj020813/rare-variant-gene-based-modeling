#!/bin/bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -o pipefail
export TMPDIR=${PROJECT_ROOT}/work/tmp TMP=${PROJECT_ROOT}/work/tmp
export TEMP=${PROJECT_ROOT}/work/tmp
mkdir -p "$TMPDIR"
FS=${PROJECT_ROOT}/work/fset
OUT=$FS/f3f6
LOG=$OUT/logs
mkdir -p "$OUT" "$LOG" "$OUT/geno"
PY=python3
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
export NUMEXPR_NUM_THREADS=4 PYTHONPATH=$OUT
cd "$OUT" || exit 1
R="nice -n 19 ionice -c3 $PY"

echo "[$(date +%F_%T)] START f3f6 chain" | tee -a "$LOG/run.log"
echo "[$(date +%F_%T)] loadavg $(cat /proc/loadavg)" | tee -a "$LOG/run.log"

step () {
  local n=$1; shift
  if [ -f "$OUT/$n.done" ]; then
    echo "[$(date +%F_%T)] SKIP $n (done)" | tee -a "$LOG/run.log"; return 0
  fi
  echo "[$(date +%F_%T)] BEGIN $n" | tee -a "$LOG/run.log"
  "$@" > "$LOG/$n.log" 2>&1
  local rc=$?
  echo "[$(date +%F_%T)] END $n rc=$rc" | tee -a "$LOG/run.log"
  if [ $rc -ne 0 ]; then echo "BLOCKED_$n rc=$rc" > "$OUT/$n.failed"; fi
  return $rc
}

step prep $R "$OUT/f3f6_prep.py" || exit 1

if [ ! -f "$OUT/geno.done" ]; then
  echo "[$(date +%F_%T)] BEGIN geno" | tee -a "$LOG/run.log"
  for c in 1 2 3 6 4 5 7 8 11 12 17 19 9 10 16 20 14 15 18 13 21 22; do \
    echo $c; done | xargs -P 8 -I{} bash -c \
    "nice -n 19 ionice -c3 $PY $OUT/f46_geno.py {} > $LOG/geno_chr{}.log 2>&1 \
     || echo BLOCKED_geno_chr{} >> $OUT/geno.failed"
  n=$(ls "$OUT"/f46_geno_chr*.done 2>/dev/null | wc -l)
  echo "[$(date +%F_%T)] END geno chr_done=$n" | tee -a "$LOG/run.log"
  [ "$n" -ge 20 ] && echo "$n" > "$OUT/geno.done"
fi

step f3 $R "$OUT/f3_coherence.py"
step f4 $R "$OUT/f4_ld.py"
step f5 $R "$OUT/f5_freq.py"
step f6 $R "$OUT/f6_usable.py"

echo "[$(date +%F_%T)] CHAIN COMPLETE" | tee -a "$LOG/run.log"
echo done > "$OUT/chain.done"
