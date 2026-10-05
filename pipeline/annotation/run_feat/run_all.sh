#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
exec 9>${PROJECT_ROOT}/work/run_feat/.feat.lock
flock -n 9 || { echo 'another feature build is running'; exit 0; }
PY=python3
OUT=${PROJECT_ROOT}/work/ref/features
one(){ N=$1; [ -e "$OUT/chr$N.done" ] && { echo "[chr$N] skip"; return 0; }; 
  $PY ${PROJECT_ROOT}/work/run_feat/build_features.py $N > ${PROJECT_ROOT}/work/run_feat/chr$N.log 2>&1 \
    && echo "[chr$N] OK" || echo "[chr$N] FAIL"; }
export -f one; export OUT PY
seq 1 22 | xargs -P 6 -I{} bash -c 'one {}'
echo "FEATURES_COMPLETE done=$(ls $OUT/*.done 2>/dev/null | wc -l)/22"
