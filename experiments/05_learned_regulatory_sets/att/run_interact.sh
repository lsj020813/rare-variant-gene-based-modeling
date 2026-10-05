#!/bin/bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -u
ATT=${PROJECT_ROOT}/work/fset/att
OUT=$ATT/interact
PY=python3
mkdir -p $OUT
cd $ATT
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
export TMPDIR=${PROJECT_ROOT}/work/tmp
export MAXV=${MAXV:-50} MINV=${MINV:-8} NPERM=${NPERM:-10}
N=${NWORK:-10}
echo "INTERACT_START MAXV=$MAXV N=$N $(date -Is)" | tee -a $OUT/run.log
for ch in $(seq 1 22); do
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 5; done
  ( nice -n 19 ionice -c3 $PY interact_ceiling.py $ch \
      > $OUT/chr${ch}.log 2>&1 ; \
    echo "chr${ch} exit=$?" >> $OUT/exit.log ) &
done
wait
echo "INTERACT_DONE $(date -Is)" | tee -a $OUT/run.log
nice -n 19 $PY aggregate_interact.py > $OUT/aggregate.log 2>&1
echo "AGGREGATE exit=$? $(date -Is)" | tee -a $OUT/run.log
touch $OUT/interact.done
