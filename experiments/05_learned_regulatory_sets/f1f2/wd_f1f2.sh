#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
W=${PROJECT_ROOT}/work/fset/f1f2
while true; do
  echo "$(date -Is) load=$(cut -d' ' -f1 /proc/loadavg) memavailG=$(free -g | awk 'NR==2{print $7}') procs=$(pgrep -cf '[f]1_tendency|[f]2_model') rss_gb=$(ps -o rss= -C python 2>/dev/null | awk '{s+=$1} END{printf "%.1f", s/1048576}') done=$(ls $W/*/*.done 2>/dev/null | wc -l) last=$(tail -1 $W/logs/run.log 2>/dev/null | cut -c12-60)" >> $W/logs/watchdog_f1f2.log
  [ -f $W/ALL.done ] && exit 0
  sleep 60
done
