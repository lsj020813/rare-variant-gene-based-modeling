#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
W=${PROJECT_ROOT}/work/fset
while true; do
  echo "$(date -Is) load=$(cut -d' ' -f1 /proc/loadavg) memfreeG=$(free -g | awk 'NR==2{print $7}') uni=$(ls $W/out/uni_chr*.tsv.gz 2>/dev/null | wc -l)/22 proc=$(pgrep -cf f_universe)" >> $W/logs/watchdog_f.log
  sleep 60
done
