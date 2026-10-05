#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -u
LIST="${1:?}"; P="${2:?}"
cd ${PROJECT_ROOT}/work/run_annot_bbj
PY=python3
run1(){ N=$1
  bash -c "ulimit -v 4000000; bash gpn_bbj.sh $N" > gpn_bbj_$N.log 2>&1; echo "[chr$N] gpn rc=$?"
  bash -c "ulimit -v 8000000; $PY t2_bbj.py $N" > t2_bbj_$N.log 2>&1; echo "[chr$N] t2 rc=$?"; }
export -f run1; export PY
printf '%s\n' $LIST | xargs -P $P -I{} bash -c 'run1 {}'
echo GT_POOL_DONE
