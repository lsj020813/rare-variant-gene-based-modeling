#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
exec 9>${PROJECT_ROOT}/work/run_b6/.b6.lock
flock -n 9 || { echo 'another b6 join running'; exit 0; }
one(){ local N=$1
  [ -e ${PROJECT_ROOT}/work/ref/b6_cards/chr$N.done ] && { echo "[chr$N] skip"; return 0; }
  ${PYTHON:-python3} \
    ${PROJECT_ROOT}/work/run_b6/b6_card_join.py chr$N \
    > ${PROJECT_ROOT}/work/run_b6/chr$N.log 2>&1 && echo "[chr$N] OK" || echo "[chr$N] FAIL"; }
export -f one
seq 1 22 | xargs -P 6 -I{} bash -c 'one "$@"' _ {}
echo "B6_ALL_COMPLETE done=$(ls ${PROJECT_ROOT}/work/ref/b6_cards/*.done 2>/dev/null | wc -l)/22"
