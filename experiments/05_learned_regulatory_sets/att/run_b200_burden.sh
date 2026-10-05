#!/bin/bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -u
ATT=${PROJECT_ROOT}/work/fset/att
OUT=$ATT/b200
PY=python3
mkdir -p $OUT/logs $OUT/private/burden
cd $ATT
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export TMPDIR=${PROJECT_ROOT}/work/tmp
N=12
echo "BURDEN_START $(date -Is)" | tee -a $OUT/logs/run.log
for ch in $(seq 1 22); do
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 5; done
  ( nice -n 19 ionice -c3 $PY att_burden_b200.py $ch \
      > $OUT/logs/burden_chr${ch}.log 2>&1 ; \
    echo "chr${ch} exit=$?" >> $OUT/logs/burden_exit.log ) &
done
wait
echo "ALL_CHR_DONE $(date -Is)" | tee -a $OUT/logs/run.log
nice -n 19 ionice -c3 $PY att_consolidate_b200.py \
  > $OUT/logs/consolidate.log 2>&1
echo "CONSOLIDATE exit=$? $(date -Is)" | tee -a $OUT/logs/run.log
touch $OUT/burden_stage.done
