#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
O=${PROJECT_ROOT}/work/ref/saige_step2_bwg
for pass in 1 2 3 4 5; do
  while ! flock -n 8 8>${PROJECT_ROOT}/work/run_bwg2/.lock; do sleep 120; done
  exec 8>&-
  n=$(ls $O/*.done 2>/dev/null | wc -l)
  echo "pass $pass: $n/88 done"
  [ "$n" -ge 88 ] && { echo SUPERVISE_ALL_DONE; exit 0; }
  export TMPDIR=${PROJECT_ROOT}/work/tmp POOL=60
  bash ${PROJECT_ROOT}/work/run_bwg2/assoc.sh >> ${PROJECT_ROOT}/work/run_bwg2/assoc_run2.log 2>&1
done
echo SUPERVISE_EXHAUSTED
