#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -u
one(){ N=$1; O=${PROJECT_ROOT}/work/ref/annot/ds/chr$N.ds.npz; [ -s $O.done ] && { echo "[chr$N] skip"; return 0; }
  cd ${PROJECT_ROOT}/work/run_annot; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
  /usr/bin/time -f "[chr$N] rss=%MKB wall=%es" python3 ds_store.py $N > ds_$N.log 2>&1 || echo "[chr$N] FAIL"; tail -1 ds_$N.log; }
export -f one
printf '%s\n' 1 2 3 6 4 7 5 12 11 10 8 9 17 16 19 14 20 18 13 15 21 | xargs -P 4 -I{} bash -c 'one {}'
echo DSPOOL_DONE
