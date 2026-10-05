#!/bin/bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -u
ATT=${PROJECT_ROOT}/work/fset/att
OUT=$ATT/oracle_hq
PY=python3
mkdir -p $OUT/logs $OUT/private/burden
cd $ATT
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
export TMPDIR=${PROJECT_ROOT}/work/tmp
export DOMSTEP=${DOMSTEP:-2} NPERM=${NPERM:-10}
N=${NWORK:-10}
echo "ORACLE_HQ_START DOMSTEP=$DOMSTEP NPERM=$NPERM N=$N $(date -Is)" | tee -a $OUT/logs/run.log
for ch in $(seq 1 22); do
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 5; done
  ( nice -n 19 ionice -c3 $PY oracle_domain_hq.py $ch \
      > $OUT/logs/chr${ch}.log 2>&1 ; \
    echo "chr${ch} exit=$?" >> $OUT/logs/exit.log ) &
done
wait
echo "ORACLE_HQ_DONE $(date -Is)" | tee -a $OUT/logs/run.log
nice -n 19 $PY aggregate_oracle_hq.py > $OUT/logs/aggregate.log 2>&1
echo "AGGREGATE exit=$? $(date -Is)" | tee -a $OUT/logs/run.log
touch $OUT/oracle.done
