#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -o pipefail
export TMPDIR=${PROJECT_ROOT}/work/tmp TMP=${PROJECT_ROOT}/work/tmp TEMP=${PROJECT_ROOT}/work/tmp
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 NUMEXPR_NUM_THREADS=4
W=${PROJECT_ROOT}/work/fset/att
OUT=$W/cvmin
PY=python3
mkdir -p $OUT/logs
ulimit -v $((32*1024*1024))
echo "$(date -Is) START att_fit_cvmin (nshuf=20)" >> $OUT/logs/run.log
cd $W && nice -n 19 ionice -c3 $PY att_fit_cvmin.py 20 > $OUT/logs/fit.log 2>&1
rc=$?
echo "$(date -Is) END att_fit_cvmin exit=$rc" >> $OUT/logs/run.log
exit $rc
