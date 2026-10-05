#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -u
LIST="${1:?chr list}"; P="${2:?pool size}"
cd ${PROJECT_ROOT}/work/run_annot_bbj
run1(){ N=$1; O=${PROJECT_ROOT}/work/ref/annot/extract_bbj/chr$N.annot.tsv
  [ -s $O.done ] && [ -s ${PROJECT_ROOT}/work/ref/annot/extract_bbj/chr$N.cadd.tsv.done ] && { echo "[chr$N] skip"; return 0; }
  ionice -c3 nice -n10 bash -c "ulimit -v 4000000; bash ext_bbj.sh $N" > ext_bbj_$N.log 2>&1; echo "[chr$N] rc=$?"; }
export -f run1
printf '%s\n' $LIST | xargs -P $P -I{} bash -c 'run1 {}'
echo POOL_BBJ_DONE
