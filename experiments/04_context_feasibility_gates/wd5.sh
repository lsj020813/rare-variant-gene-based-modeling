#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
W=${PROJECT_ROOT}/work/gate1
while true; do
  echo "$(date -Is) load=$(cut -d' ' -f1-3 /proc/loadavg) memfreeG=$(free -g | awk 'NR==2{print $7}') files=$(ls $W/out/universe_chr*.tsv 2>/dev/null | wc -l) proc=$(pgrep -cf 'gate1_universe_census')" >> $W/logs/watchdog.log
  sleep 60
done
