#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
W=${PROJECT_ROOT}/work/gate2
while true; do
  echo "$(date -Is) load=$(cut -d' ' -f1 /proc/loadavg) memfreeG=$(free -g | awk 'NR==2{print $7}') g2=$(ls $W/out/g2_chr*.json 2>/dev/null | wc -l)/22 proc=$(pgrep -cf gate2_context_v2)" >> $W/logs/watchdog_g2.log
  sleep 60
done
