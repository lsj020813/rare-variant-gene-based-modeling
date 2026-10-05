#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -o pipefail
LIST=$1; P=${2:-2}; ARMS=${3:-gate,phi,perm}; NPERM=${4:-20}; TAG=${5:-}
W=${PROJECT_ROOT}/work/run_trackA; cd $W
LOG=$W/logs/pool${TAG:+_$TAG}.log; log(){ echo "$(date '+%F %T') $*" >> $LOG; }
log "POOL START list=$LIST P=$P arms=$ARMS nperm=$NPERM n=$(wc -l < $LIST)"
cnt(){ if [ -n "$TAG" ]; then ls $W/fm/*.$TAG.done 2>/dev/null | wc -l; else ls $W/fm/*.done 2>/dev/null | grep -v -E '\.(gate|smoke)\.done$' | grep -v POOL | wc -l; fi; }
export -f cnt
run1(){ T=$1; RID=$2; bash $W/scripts/run_region.sh $T $RID "$ARMS" "$NPERM" "$TAG" >> $W/logs/pool_regions${TAG:+_$TAG}.log 2>&1; echo "$(date '+%F %T') $RID rc=$?" >> $W/logs/pool${TAG:+_$TAG}.log
  ND=$(cnt)
  if [ $(( ND % 20 )) -eq 0 ] && [ -z "$TAG" ]; then echo "$(date '+%F %T') [fm] 완료 구역 $ND: $(bash $W/scripts/partial_summary.sh 2>/dev/null)" >> $W/PROGRESS.md; fi; }
export -f run1; export W ARMS NPERM TAG
cut -f1,2 $LIST | xargs -P $P -L 1 bash -c 'run1 $0 $1'
ND=$(cnt); NS=$(ls $W/fm/*.skip 2>/dev/null | wc -l)
log "POOL END done=$ND skip=$NS"; echo "$(date '+%F %T') [fm] 풀 종료 done=$ND skip=$NS: $(bash $W/scripts/partial_summary.sh 2>/dev/null)" >> $W/PROGRESS.md
touch $W/fm/POOL${TAG:+.$TAG}.done
