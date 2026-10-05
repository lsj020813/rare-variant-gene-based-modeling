#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -u
SC=$1; P=$2; UL=$3; TAG=$4; cd ${PROJECT_ROOT}/work/run_annot
export TMPDIR=${PROJECT_ROOT}/work/tmp
echo "[$TAG] $(date '+%m-%d %H:%M') start pool=$P ulimit=${UL}kb"
printf '%s\n' 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 | \
  xargs -P $P -I{} bash -c "ulimit -v $UL; echo 800 > /proc/self/oom_score_adj 2>/dev/null; bash $SC {} >> ${TAG}.chr{}.log 2>&1; echo \"[$TAG] chr{} exit \$?\""
echo "[$TAG] $(date '+%m-%d %H:%M') POOL_DONE done=$(ls ${PROJECT_ROOT}/work/ref/annot/extract/*.cadd.tsv.done ${PROJECT_ROOT}/work/ref/common05/*.done 2>/dev/null | wc -l)"
