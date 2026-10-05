#!/bin/bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -u
ATT=${PROJECT_ROOT}/work/fset/att
PY=python3
mkdir -p $ATT/logs $ATT/private/burden
cd $ATT
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export TMPDIR=${PROJECT_ROOT}/work/tmp
N=12
for ch in $(seq 1 22); do
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 5; done
  ( nice -n 19 ionice -c3 $PY att_burden.py $ch \
      > $ATT/logs/burden_chr${ch}.log 2>&1 ; \
    echo "chr${ch} exit=$?" >> $ATT/logs/burden_exit.log ) &
done
wait
echo "ALL_CHR_DONE $(date -Is)" | tee -a $ATT/logs/run.log
nice -n 19 ionice -c3 $PY att_consolidate.py \
  > $ATT/logs/consolidate.log 2>&1
echo "CONSOLIDATE exit=$? $(date -Is)" | tee -a $ATT/logs/run.log
touch $ATT/burden_stage.done
