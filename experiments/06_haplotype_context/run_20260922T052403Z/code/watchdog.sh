#!/bin/bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
cd ${PROJECT_ROOT}/work/haplotype_context_20260922_v1/run_20260922T052403Z
LOG=logs/watchdog.log
while true; do
  TS=$(date -u +'%Y-%m-%dT%H:%M:%SZ')
  LOAD=$(cut -d' ' -f1-3 /proc/loadavg)
  MEM=$(free -g | awk 'NR==2{print "used="$3"G avail="$7"G"}')
  NPROC=$(pgrep -u "$USER" -f '[h]cf.build_states' | wc -l)
  NTILES=$(grep -ch . logs/../states/worker*.log 2>/dev/null | paste -sd+ | bc 2>/dev/null)
  echo "$TS load=$LOAD $MEM my_workers=$NPROC worker_log_lines=${NTILES:-0}" >> $LOG
  if [ -f logs/build_states.done ]; then
    echo "$TS DONE marker present rc=$(cat logs/build_states.done)" >> $LOG; break
  fi
  if [ "$NPROC" -eq 0 ]; then
    if [ ! -f logs/build_states.done ]; then
      sleep 20
      NPROC2=$(pgrep -u "$USER" -f '[h]cf.build_states' | wc -l)
      if [ "$NPROC2" -eq 0 ] && [ ! -f logs/build_states.done ]; then
        echo "$TS WATCHDOG_NO_WORKERS_NO_DONE_MARKER" >> $LOG; break
      fi
    fi
  fi
  sleep 60
done
