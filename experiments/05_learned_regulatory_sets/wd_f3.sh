#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
W=${PROJECT_ROOT}/work/fset; O=$W/primary
while true; do
  echo "$(date -Is) load=$(cut -d' ' -f1 /proc/loadavg) memavailG=$(free -g | awk 'NR==2{print $7}') feat=$(ls $O/feat_chr*.tsv.gz 2>/dev/null | wc -l)/22 vset=$(ls $O/vset_chr*.tsv.gz 2>/dev/null | wc -l)/22 gpn=$([ -f $O/gpn.done ] && echo done || echo run) featproc=$(pgrep -cf '[f]_features_v3') gpnproc=$(pgrep -cf '[f]_gpn_v2|[t]abix -R') cset=$([ -f $O/cset.done ] && echo done || echo -)" >> $W/logs/watchdog_f3.log
  sleep 60
done
