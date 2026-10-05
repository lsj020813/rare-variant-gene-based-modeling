#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -u
run1(){ N=$1; O=${PROJECT_ROOT}/work/ref/annot/gpn/chr$N.gpn.tsv
  [ -s $O.done ] && { echo "[chr$N] skip"; return 0; }
  bash -c "ulimit -v 4000000; bash gpn1.sh $N" > gpn$N.log 2>&1; }
export -f run1
printf '%s\n' 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 | xargs -P 3 -I{} bash -c 'run1 {}'
echo GPOOL_DONE
