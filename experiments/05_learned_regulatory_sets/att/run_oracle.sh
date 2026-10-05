#!/bin/bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -u
ATT=${PROJECT_ROOT}/work/fset/att
OUT=$ATT/oracle
PY=python3
mkdir -p $OUT/logs $OUT/private/burden
cd $ATT
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
export TMPDIR=${PROJECT_ROOT}/work/tmp
export DOMSTEP=${DOMSTEP:-1} NPERM=${NPERM:-10}
N=12
echo "ORACLE_START DOMSTEP=$DOMSTEP NPERM=$NPERM $(date -Is)" | tee -a $OUT/logs/run.log
for ch in $(seq 1 22); do
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 5; done
  ( nice -n 19 ionice -c3 $PY oracle_domain.py $ch \
      > $OUT/logs/chr${ch}.log 2>&1 ; \
    echo "chr${ch} exit=$?" >> $OUT/logs/exit.log ) &
done
wait
echo "ORACLE_DONE $(date -Is)" | tee -a $OUT/logs/run.log
nice -n 19 $PY aggregate_oracle.py > $OUT/logs/aggregate.log 2>&1
echo "AGGREGATE exit=$? $(date -Is)" | tee -a $OUT/logs/run.log
touch $OUT/oracle.done
