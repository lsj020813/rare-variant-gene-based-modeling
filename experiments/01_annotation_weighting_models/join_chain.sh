#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -u; cd ${PROJECT_ROOT}/work/run_band15; export OMP_NUM_THREADS=2 TMPDIR=${PROJECT_ROOT}/work/tmp
R=${PROJECT_ROOT}/work/ref15; PY=python3
echo "[chain] $(date '+%m-%d %H:%M') stage1 groupfiles"
gf(){ N=$1; [ -s $R/groupfiles_bwg/chr$N.summary.json ] && return 0; $PY bwg.py $N > bwg_$N.log 2>&1 || echo "[chr$N] BWG FAIL"; }
export -f gf; export R PY; printf '%s\n' $(seq 1 22) | xargs -P 4 -I{} bash -c 'gf {}'
n=$(ls $R/groupfiles_bwg/*.summary.json | wc -l); [ "$n" -eq 22 ] || { echo "[chain] GATE FAIL groupfiles $n/22"; exit 2; }
echo "[chain] $(date '+%H:%M') stage2 map"
[ -s $R/annot/map38/chr1.map.tsv ] && [ $(ls $R/annot/map38/*.map.tsv | wc -l) -eq 22 ] || $PY map.py > map_full.log 2>&1
n=$(ls $R/annot/map38/*.map.tsv | wc -l); [ "$n" -eq 22 ] || { echo "[chain] GATE FAIL map $n/22"; exit 3; }
echo "[chain] $(date '+%H:%M') stage3 favor (pool 6, ionice)"
ex(){ N=$1; [ -s $R/annot/extract/chr$N.annot.tsv.done ] && return 0; ionice -c3 nice -n10 bash -c "ulimit -v 4000000; bash ext.sh $N" > ext$N.log 2>&1 || echo "[chr$N] EXT FAIL"; }
export -f ex; printf '%s\n' 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 | xargs -P 6 -I{} bash -c 'ex {}'
n=$(ls $R/annot/extract/*.done | wc -l); [ "$n" -eq 22 ] || { echo "[chain] GATE FAIL favor $n/22"; exit 4; }
echo "[chain] $(date '+%H:%M') stage4 gpn (pool 3)"
gp(){ N=$1; [ -s $R/annot/gpn/chr$N.gpn.tsv.done ] && return 0; bash -c "ulimit -v 4000000; bash gpn1.sh $N" > gpn$N.log 2>&1 || echo "[chr$N] GPN FAIL"; }
export -f gp; printf '%s\n' $(seq 1 22) | xargs -P 3 -I{} bash -c 'gp {}'
n=$(ls $R/annot/gpn/*.done | wc -l); [ "$n" -eq 22 ] || { echo "[chain] GATE FAIL gpn $n/22"; exit 5; }
echo "[chain] $(date '+%H:%M') stage5 t2/t4/assemble (pool 2)"
fm(){ N=$1; [ -s $R/annot/t2/chr$N.t2.tsv.done ] || bash -c "ulimit -v 8000000; $PY t2_link.py $N" > t2_$N.log 2>&1 || { echo "[chr$N] T2 FAIL"; return 1; }
  bash t4one.sh $N > t4_$N.log 2>&1 || { echo "[chr$N] T4 FAIL"; return 1; }
  [ -s $R/annot/fm/chr$N.fm.tsv ] || bash -c "ulimit -v 4000000; $PY assemble.py $N" > fm_$N.log 2>&1 || { echo "[chr$N] FM FAIL"; return 1; }; echo "[chr$N] fm ok"; }
export -f fm; printf '%s\n' $(seq 1 22) | xargs -P 2 -I{} bash -c 'fm {}'
n=$(ls $R/annot/fm/*.fm.tsv | wc -l); echo "[chain] $(date '+%H:%M') fm files $n/22"
[ "$n" -eq 22 ] && echo BAND15_JOIN_CHAIN_OK || echo "[chain] GATE FAIL fm $n/22"
