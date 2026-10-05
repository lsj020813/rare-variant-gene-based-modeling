#!/usr/bin/env bash
PAT="$1"; LOGF="$2"; LOADMAX=${3:-40}; MEMMIN=${4:-20}
echo "START $(date -Is) pat=$PAT loadmax=$LOADMAX memmin=${MEMMIN}G" >> "$LOGF"
n=0; stopped=0; miss=0
while true; do
  n=$((n+1))
  pids=$(pgrep -f "$PAT" | tr '\n' ' ')
  if [ -z "$pids" ]; then
    miss=$((miss+1))
    if [ "$miss" -ge 3 ]; then echo "$(date -Is) n=$n target gone -- watchdog exit" >> "$LOGF"; exit 0; fi
  else
    miss=0
  fi
  l1=$(cut -d' ' -f1 /proc/loadavg)
  memav=$(awk '/MemAvailable/{printf "%d", $2/1048576}' /proc/meminfo)
  rss=$(ps -o rss= -p $(echo $pids | tr ' ' ',') 2>/dev/null | awk '{s+=$1} END{printf "%.1f", s/1048576}')
  over=$(awk -v a="$l1" -v b="$LOADMAX" 'BEGIN{print (a>b)?1:0}')
  low=$(awk -v a="$memav" -v b="$MEMMIN" 'BEGIN{print (a<b)?1:0}')
  if [ "$over" = "1" ] || [ "$low" = "1" ]; then
    if [ "$stopped" = "0" ] && [ -n "$pids" ]; then kill -STOP $pids 2>/dev/null; stopped=1
      echo "$(date -Is) n=$n STOP load1=$l1 memavail=${memav}G rss=${rss}G" >> "$LOGF"; fi
  else
    if [ "$stopped" = "1" ] && [ -n "$pids" ]; then kill -CONT $pids 2>/dev/null; stopped=0
      echo "$(date -Is) n=$n CONT load1=$l1 memavail=${memav}G rss=${rss}G" >> "$LOGF"; fi
  fi
  echo "$(date -Is) n=$n load1=$l1 memavail=${memav}G rss=${rss}G stopped=$stopped" >> "$LOGF"
  sleep 60
done
