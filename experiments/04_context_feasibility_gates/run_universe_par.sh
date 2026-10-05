#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
W=${PROJECT_ROOT}/work/gate1; mkdir -p $W/out $W/logs
echo "START $(date -Is) pid=$$"
JOBS=$W/logs/joblist.txt; : > $JOBS
for N in 1 2 4 3 6 5 8 7 12 11 10 9 17 16 19 14 18 20 13 22 21 15; do
  echo "$N ${PROJECT_ROOT}/work/ref/orig_index/chr$N.vcf.gz orig" >> $JOBS
  echo "$N ${PROJECT_ROOT}/work/ref/common05/chr$N.maf05.vcf.gz common05" >> $JOBS
done
cat $JOBS | xargs -P 8 -L 1 bash $W/gate1_universe_par.sh
echo "END $(date -Is)"
ls $W/out/universe_orig_chr*.tsv | wc -l > $W/out/universe_orig.count
echo "universe census complete $(date -Is) orig=$(ls $W/out/universe_orig_chr*.tsv 2>/dev/null | wc -l) common05=$(ls $W/out/universe_common05_chr*.tsv 2>/dev/null | wc -l)" > $W/out/universe.done
