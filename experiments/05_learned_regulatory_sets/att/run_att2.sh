#!/bin/bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -u
ATT=${PROJECT_ROOT}/work/fset/att
PY=python3
cd $ATT
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export TMPDIR=${PROJECT_ROOT}/work/tmp
N=12
ORDER="1 11 5 2 7 19 12 6 16 17 20 3 9 14 10 4 22 21 8 18 15 13"
for ch in $ORDER; do
  if [ -f $ATT/burden_chr${ch}.done ]; then continue; fi
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 10; done
  ( nice -n 19 ionice -c3 $PY att_burden.py $ch > $ATT/logs/burden_chr${ch}.log 2>&1
    echo "chr${ch} exit=$? $(date -Is)" >> $ATT/logs/burden_exit.log ) &
  sleep 2
done
wait
:
echo "ALL_CHR_DONE $(date -Is)" >> $ATT/logs/run.log
nice -n 19 ionice -c3 $PY att_consolidate.py > $ATT/logs/consolidate.log 2>&1
echo "CONSOLIDATE exit=$? $(date -Is)" >> $ATT/logs/run.log
nice -n 19 ionice -c3 $PY att_bbj.py > $ATT/logs/bbj.log 2>&1
echo "BBJ exit=$? $(date -Is)" >> $ATT/logs/run.log
touch $ATT/burden_stage.done
