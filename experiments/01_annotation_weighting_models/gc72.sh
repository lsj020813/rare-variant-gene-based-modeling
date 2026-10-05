#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
cd ${PROJECT_ROOT}/work/run_band15
export OMP_NUM_THREADS=4 CUDA_VISIBLE_DEVICES=
PY=python3
for arm in cluster spline linear; do
  echo "=== $arm ==="
  (ulimit -v 40000000; $PY l1_train_v72.py --smoke 22 --arm $arm --gradcheck-only 2>&1 | grep -E 'design|cluster K|GRADCHECK|PD ok|GATE|Error|Traceback|S2 indicators|DONE|multiple chrom|mismatch' | cut -c1-260)
done
echo GC72_ALL_DONE
