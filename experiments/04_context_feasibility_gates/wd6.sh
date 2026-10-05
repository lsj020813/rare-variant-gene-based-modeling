#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
W=${PROJECT_ROOT}/work/gate1
while true; do
  L=$(cut -d' ' -f1 /proc/loadavg); M=$(free -g | awk 'NR==2{print $7}')
  U=$(ls $W/out/universe_*_chr*.tsv 2>/dev/null | wc -l); G=$(ls $W/out/window_summary_chr*.tsv 2>/dev/null | wc -l)
  C=$(ls $W/out/ctx_chr*.json 2>/dev/null | wc -l)
  P=$(pgrep -cf 'gate1_universe_par|gate1_windows|gate1_context' 2>/dev/null)
  echo "$(date -Is) load=$L memfreeG=$M universe=$U/44 windows=$G/22 ctx=$C proc=$P" >> $W/logs/watchdog6.log
  awk -v l="$L" 'BEGIN{ if (l+0 > 60) print "ALERT high load" }' >> $W/logs/watchdog6.log
  sleep 60
done
