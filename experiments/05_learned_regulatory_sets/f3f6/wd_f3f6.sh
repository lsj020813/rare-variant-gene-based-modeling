#!/bin/bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
PAT=${PAT:-'[f]set/f3f6/f'}
LOG=${PROJECT_ROOT}/work/fset/f3f6/logs/watchdog_f3f6.log
NCORE=$(grep -c ^processor /proc/cpuinfo)
LIMIT=${LIMIT:-40}
echo "[$(date +%F_%T)] WD START ncore=$NCORE limit=$LIMIT" >> "$LOG"
for i in $(seq 1 60); do
  pgrep -f "$PAT" > /dev/null && break
  sleep 5
done
if ! pgrep -f "$PAT" > /dev/null; then
  echo "[$(date +%F_%T)] WD WARN target never appeared" >> "$LOG"; exit 2
fi
stopped=0; n=0
while true; do
  pids=$(pgrep -f "$PAT")
  [ -z "$pids" ] && { echo "[$(date +%F_%T)] WD target gone n=$n -- exit" \
    >> "$LOG"; break; }
  n=$((n+1))
  l1=$(awk '{print $1}' /proc/loadavg)
  free_g=$(free -g | awk '/^Mem:/{print $7}')
  ours=$(ps -o rss= -p $(echo $pids | tr ' ' ',') 2>/dev/null \
    | awk '{s+=$1} END{print int(s/1048576)}')
  over=$(awk -v a="$l1" -v b="$LIMIT" 'BEGIN{print (a>b)?1:0}')
  low=$(awk -v a="$free_g" 'BEGIN{print (a<20)?1:0}')
  if [ "$over" = 1 ] || [ "$low" = 1 ]; then
    if [ "$stopped" = 0 ]; then
      kill -STOP $pids 2>/dev/null; stopped=1
      echo "[$(date +%F_%T)] WD STOP n=$n load1=$l1 availGB=$free_g" >> "$LOG"
    fi
  else
    if [ "$stopped" = 1 ]; then
      kill -CONT $pids 2>/dev/null; stopped=0
      echo "[$(date +%F_%T)] WD CONT n=$n load1=$l1 availGB=$free_g" >> "$LOG"
    fi
  fi
  [ $((n % 12)) = 1 ] && echo \
    "[$(date +%F_%T)] WD n=$n load1=$l1 availGB=$free_g ourRSSg=$ours" \
    >> "$LOG"
  sleep 30
done
