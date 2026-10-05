#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -o pipefail
W=${PROJECT_ROOT}/work/run_trackA; A=${PROJECT_ROOT}/work/ref/annot; S=$W/scripts
PY=python3
export TMPDIR=${PROJECT_ROOT}/work/tmp TMP=${PROJECT_ROOT}/work/tmp TEMP=${PROJECT_ROOT}/work/tmp
cd $W; mkdir -p logs annot_trA
LOG=$W/logs/annot_trA.log; log(){ echo "$(date '+%F %T') $*" >> $LOG; }; prog(){ echo "$(date '+%F %T') $*" >> $W/PROGRESS.md; }
log "ANNOT_TRA START"
if [ ! -s annot_trA/missing_keys37.txt.done ]; then $PY $S/missing_keys.py > logs/missing_keys.log 2>&1 || { log "GATE FAIL: missing_keys"; exit 2; }; fi
log "missing_keys: $(tail -1 logs/missing_keys.log)"
if [ ! -s $A/map38_trA/MAP_DONE ]; then ( ulimit -v 16000000; $PY $S/map_trA.py > logs/map_trA.log 2>&1 ) || { log "GATE FAIL: map_trA"; exit 3; }; fi
log "map_trA: $(grep -o '{.*}' logs/map_trA.log | tail -1 | cut -c1-300)"
CHRS=""; SKIPC=""
for f in $(ls $A/map38_trA/chr*.map.tsv.done | sort -V); do N=$(echo $f | sed 's/.*chr\([0-9]*\)\.map.*/\1/'); C=$(cut -d' ' -f2 $f); if [ "$C" -ge 20 ]; then CHRS="$CHRS $N"; else SKIPC="$SKIPC chr$N:$C"; fi; done
CHRS=$(echo $CHRS | tr ' ' '\n' | sort -n)
log "chrs: $(echo $CHRS | tr '\n' ' ') | skipped small buckets:$SKIPC"
prog "[annot] missing_keys·map 완료: $(tail -1 logs/missing_keys.log) / map $(grep -o '"hit": [0-9]*' logs/map_trA.log | tail -1)"
run_ext(){ N=$1; O=$A/extract_trA/chr$N.annot.tsv
  [ -s $O.done ] && [ -s $A/extract_trA/chr$N.cadd.tsv.done ] && { echo "[chr$N] skip"; return 0; }
  ionice -c3 nice -n10 bash -c "ulimit -v 4000000; bash $S/ext_trA.sh $N" > $W/logs/ext_trA_$N.log 2>&1; echo "[chr$N] rc=$?"; }
export -f run_ext; export A S W
( printf '%s\n' $CHRS | sort -n | xargs -P 4 -I{} bash -c 'run_ext {}' > logs/ext_pool.log 2>&1; echo EXT_POOL_DONE >> logs/ext_pool.log ) &
for N in $CHRS; do
  [ -s $A/gpn_trA/chr$N.gpn.tsv.done ] || bash $S/gpn_trA.sh $N > logs/gpn_trA_$N.log 2>&1 || log "GATE FAIL: gpn chr$N"
  [ -s $A/t2_trA/chr$N.t2.tsv.done ] || $PY $S/t2_trA.py $N > logs/t2_trA_$N.log 2>&1 || log "GATE FAIL: t2 chr$N"
done
log "gpn/t2 loop done: gpn $(ls $A/gpn_trA/*.done 2>/dev/null | wc -l) t2 $(ls $A/t2_trA/*.done 2>/dev/null | wc -l)"
wait
NC=$(echo $CHRS | wc -w); NE=$(ls $A/extract_trA/chr*.annot.tsv.done 2>/dev/null | wc -l)
log "ext pool done: $NE/$NC"; prog "[annot] FAVOR 추출 $NE/$NC 염색체 완료 ($(grep -h 'join' logs/ext_trA_*.log | sed 's/.*join //' | tr '\n' ' '))"
[ "$NE" -eq "$NC" ] || { log "GATE FAIL: ext incomplete"; exit 4; }
( ulimit -v 40000000; $PY $S/assemble_trA.py --chrs $(echo $CHRS | tr ' ' ',') > logs/assemble_trA.log 2>&1 ) || { log "GATE FAIL: assemble"; exit 5; }
log "assemble: $(grep -o '{.*}' logs/assemble_trA.log | tail -1 | cut -c1-400)"
prog "[annot] 조립 완료: $(grep -o '"rows": [0-9]*, "dropped_chr2122": [0-9]*, "dropped_bad_maf": [0-9]*, "chr38_ne_chr37": [0-9]*, "in_our_band": [0-9]*, "t1_na": [0-9]*' logs/assemble_trA.log | tail -1)"
touch annot_trA/ANNOT_TRA.done; log "ANNOT_TRA END"
