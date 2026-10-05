#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
cd ${PROJECT_ROOT}/work/run_annot
for i in $(seq 1 840); do
  for N in $(seq 1 22); do
    O=${PROJECT_ROOT}/work/ref/common05/chr$N.maf05.vcf.gz
    [ -s "$O.done" ] && continue
    [ -s "$O.tmp" ] || continue
    grep -q 'FORMAT != DS' common.chr$N.log 2>/dev/null || continue
    pgrep -f "ext_common.sh $N\$" >/dev/null && continue
    echo "[finish] $(date '+%m-%d %H:%M') chr$N"; bash ext_finish.sh $N
  done
  [ "$(ls ${PROJECT_ROOT}/work/ref/common05/*.done 2>/dev/null | wc -l)" -ge 22 ] && { echo "[finish] ALL 22 DONE"; exit 0; }
  sleep 60
done
