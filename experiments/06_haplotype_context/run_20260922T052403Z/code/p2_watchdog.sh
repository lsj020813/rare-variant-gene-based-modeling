#!/bin/bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
while [ ! -f ${PROJECT_ROOT}/work/haplotype_context_20260922_v1/run_20260922T052403Z/logs/p2_all.done ]; do
  d=$(date -u +%FT%TZ); l=$(cut -d' ' -f1-3 /proc/loadavg)
  mu=$(free -g | awk 'NR==2{print $3}'); ma=$(free -g | awk 'NR==2{print $7}')
  rss=$(ps -o rss= -C python 2>/dev/null | awk '{s+=$1} END {print int(s/1048576)}')
  t=$(ls ${PROJECT_ROOT}/work/haplotype_context_20260922_v1/run_20260922T052403Z/models/p2/bc_*.csv 2>/dev/null | wc -l); f=$(ls ${PROJECT_ROOT}/work/haplotype_context_20260922_v1/run_20260922T052403Z/models/p2/fail_*.json 2>/dev/null | wc -l)
  echo "$d load=$l mem_used_gb=$mu avail_gb=$ma python_rss_gb=$rss tiles=$t fail=$f" >> ${PROJECT_ROOT}/work/haplotype_context_20260922_v1/run_20260922T052403Z/logs/p2_watchdog.log
  sleep 60
done
