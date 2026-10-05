#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
K=$1; shift; CHRS="$*"
W=${PROJECT_ROOT}/work/run_ourfm; G=${PROJECT_ROOT}/work/ref/gwas05
LOG=$W/logs/pool.log
DEADLINE=$(date -d "${FM_DEADLINE:-+48 hours}" +%s)
log(){ echo "$(date '+%F %T') $*" >> $LOG; }
log "POOL START K=$K chrs=$CHRS deadline=$(date -d @$DEADLINE '+%F %T')"
declare -A started
while true; do
  now=$(date +%s); remaining=0
  for N in $CHRS; do
    [ -s $W/logs/chr$N.done ] && continue
    remaining=1
    [ "${started[$N]:-0}" -eq 1 ] && continue
    ok=1; for T in tchl htn dm lip; do [ -s $G/$T.chr$N.txt.done ] || ok=0; done
    [ $ok -eq 1 ] || continue
    [ $now -lt $DEADLINE ] || continue
    pgrep -f "scripts/fm_chr[.]sh $N\$" >/dev/null && { started[$N]=1; continue; }
    running=$(ps -o args= -p $(pgrep -d, -f 'scripts/fm_chr[.]sh') 2>/dev/null | sort -u | wc -l)
    [ "$running" -lt "$K" ] || break
    AV=$(free -g | awk 'NR==2{print $7}'); [ "$AV" -ge 100 ] || { log "램 가용 ${AV}G <100 — 대기"; break; }
    nohup bash $W/scripts/fm_chr.sh $N > /dev/null 2>&1 &
    started[$N]=1; log "LAUNCH chr$N (running=$((running+1)))"
    sleep 5
  done
  [ $remaining -eq 0 ] && { log "POOL ALL DONE"; break; }
  if [ $now -ge $DEADLINE ] && ! pgrep -f 'scripts/fm_chr[.]sh' >/dev/null; then
    log "POOL DEADLINE — 미완: $(for N in $CHRS; do [ -s $W/logs/chr$N.done ] || echo -n "$N "; done)"; break
  fi
  sleep 60
done
