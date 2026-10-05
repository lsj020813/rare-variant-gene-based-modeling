#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -u; P=$1; cd ${PROJECT_ROOT}/work/run_annot; export TMPDIR=${PROJECT_ROOT}/work/tmp
echo "[gwas] $(date '+%m-%d %H:%M') start pool=$P"
printf '%s\n' 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 | \
  xargs -P $P -I{} bash -c "ulimit -v 12000000; echo 800 > /proc/self/oom_score_adj 2>/dev/null; bash gwas_common.sh {} >> gwas.chr{}.log 2>&1; echo \"[gwas] chr{} exit \$?\""
echo "[gwas] $(date '+%m-%d %H:%M') GWAS_POOL_DONE done=$(ls ${PROJECT_ROOT}/work/ref/gwas05/*.done 2>/dev/null | wc -l)/88"
