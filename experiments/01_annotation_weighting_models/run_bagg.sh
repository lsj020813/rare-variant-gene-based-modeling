#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
exec 200>${PROJECT_ROOT}/work/run_l1/.bagg.lock; flock -n 200 || { echo LOCKED; exit 0; }
PY=python3
for N in $(seq 22 -1 1); do echo $N; done | xargs -P 8 -I{} $PY ${PROJECT_ROOT}/work/run_l1/bagg.py {} >> ${PROJECT_ROOT}/work/run_l1/bagg_all.log 2>&1
echo BAGG_ALL_DONE >> ${PROJECT_ROOT}/work/run_l1/bagg_all.log
