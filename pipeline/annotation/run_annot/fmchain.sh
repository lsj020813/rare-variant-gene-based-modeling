#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -u
one(){ N=$1; cd ${PROJECT_ROOT}/work/run_annot
  [ -s ${PROJECT_ROOT}/work/ref/annot/t2/chr$N.t2.tsv.done ] || bash -c "ulimit -v 8000000; python3 t2_link.py $N" > t2_$N.log 2>&1 || { echo "[chr$N] T2 FAIL"; return 1; }
  bash -c "ulimit -v 4000000; python3 assemble.py $N" > fm_$N.log 2>&1 || { echo "[chr$N] FM FAIL"; return 1; }
  echo "[chr$N] ok"; }
export -f one
printf '%s\n' 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 | xargs -P 2 -I{} bash -c 'one {}'
echo FMCHAIN_DONE
