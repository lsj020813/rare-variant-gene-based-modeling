#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
exec 201>${PROJECT_ROOT}/work/run_l1/.bagg2.lock; flock -n 201 || { echo LOCKED; exit 0; }
for N in $(seq 14 -1 1); do echo $N; done | xargs -P 6 -I{} bash ${PROJECT_ROOT}/work/run_l1/bagg_claim.sh {} >> ${PROJECT_ROOT}/work/run_l1/bagg_all.log 2>&1
echo BAGG_POOL2_DONE >> ${PROJECT_ROOT}/work/run_l1/bagg_all.log
