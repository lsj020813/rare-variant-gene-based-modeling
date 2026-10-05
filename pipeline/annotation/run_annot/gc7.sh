#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
cd ${PROJECT_ROOT}/work/run_annot
export OMP_NUM_THREADS=6
PY=python3
for arm in spline linear nn; do
  echo "=== $arm ==="
  $PY l1_train_v7.py --smoke 22 --arm $arm --gradcheck-only 2>&1 | grep -E 'design|GRADCHECK|PD ok|GATE|Error|Traceback|frozen|S2 indicators|nn params|DONE' | cut -c1-260
done
echo GC7_ALL_DONE
