#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
LIST=$1; P=${2:-2}; NPERM=${3:-1000}; W=${PROJECT_ROOT}/work/run_trackA; cd $W
log(){ echo "$(date '+%F %T') $*" >> $W/logs/pool_rw.log; }
log "RW POOL START list=$LIST P=$P nperm=$NPERM"
run1(){ bash $W/scripts/run_rw.sh $1 $2 $NPERM >> $W/logs/pool_rw_regions.log 2>&1; echo "$(date '+%F %T') $2 rc=$?" >> $W/logs/pool_rw.log
  ND=$(ls $W/rw/*.done 2>/dev/null | wc -l); if [ $(( ND % 20 )) -eq 0 ]; then echo "$(date '+%F %T') [rw] 재가중 완료 구역 $ND" >> $W/PROGRESS.md; fi; }
export -f run1; export W NPERM
cut -f1,2 $LIST | xargs -P $P -L 1 bash -c 'run1 $0 $1'
ND=$(ls $W/rw/*.done 2>/dev/null | wc -l); log "RW POOL END done=$ND"; echo "$(date '+%F %T') [rw] 재가중 풀 종료 done=$ND/81" >> $W/PROGRESS.md; touch $W/rw/POOL.done
