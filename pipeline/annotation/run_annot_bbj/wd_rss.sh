#!/usr/bin/env bash
PATTERN="${1:?pattern}"; LOG="${2:?log}"; BUDGET="${3:-38}"; INT="${4:-6}"
streak=0; peak=0; n=0; missing=0
echo "$(date '+%F %T') wd_rss start budget=${BUDGET}G interval=${INT}s pattern=$PATTERN" >> "$LOG"
while true; do
  pids=$(pgrep -f "$PATTERN" | tr '\n' ' ')
  if [ -z "$pids" ]; then
    missing=$((missing+1))
    if [ "$missing" -ge 5 ]; then echo "$(date '+%F %T') target gone after n=$n samples — peak=${peak}G — exit clean" >> "$LOG"; exit 0; fi
    sleep "$INT"; continue
  fi
  missing=0
  rss_kb=$(ps -o rss= -p "$(echo $pids | tr ' ' ',')" 2>/dev/null | awk '{s+=$1} END{print s+0}')
  rss_g=$(awk -v k="$rss_kb" 'BEGIN{printf "%.1f", k/1048576}')
  avail=$(awk '/MemAvailable/{printf "%d", $2/1048576}' /proc/meminfo)
  n=$((n+1)); peak=$(awk -v a="$peak" -v b="$rss_g" 'BEGIN{print (b>a)?b:a}')
  if [ "$(awk -v r="$rss_g" -v b="$BUDGET" 'BEGIN{print (r>b)?1:0}')" -eq 1 ]; then streak=$((streak+1)); else streak=0; fi
  if [ $((n % 10)) -eq 1 ] || [ "$streak" -gt 0 ]; then echo "$(date '+%F %T') n=$n rss=${rss_g}G peak=${peak}G avail=${avail}G streak=$streak" >> "$LOG"; fi
  if [ "$streak" -ge 2 ]; then echo "$(date '+%F %T') *** KILL: actual rss=${rss_g}G > ${BUDGET}G on 2 consecutive samples; pids=$pids" >> "$LOG"; kill -9 $pids; exit 1; fi
  if [ "$avail" -lt 40 ]; then echo "$(date '+%F %T') *** KILL: server MemAvailable=${avail}G < 40G; pids=$pids" >> "$LOG"; kill -9 $pids; exit 1; fi
  sleep "$INT"
done
