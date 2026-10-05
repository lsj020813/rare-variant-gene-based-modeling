#!/bin/bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
ATT=${PROJECT_ROOT}/work/fset/att
LOG=$ATT/logs/night_monitor.log
mkdir -p $ATT/logs
echo "$(date -Is) NIGHT_MONITOR START" >> $LOG
while true; do
  [ -f $ATT/night.stop ] && { echo "$(date -Is) stop marker -> exit" >> $LOG; exit 0; }
  L=$(cut -d' ' -f1-3 /proc/loadavg)
  AV=$(free -g | awk '/^Mem:/{print $7}')
  SW=$(free -m | awk '/^Swap:/{printf "%d/%d", $3, $2}')
  DF=$(df -BG --output=avail /data | tail -1 | tr -dc '0-9')
  WA=$(vmstat 1 2 | tail -1 | awk '{print $16}')
  RSS=$(ps -u user -o rss=,args= | grep -E 'att_fit|att_burden_b200|pos_control' \
        | grep -v grep | awk '{s+=$1} END{printf "%.1f", s/1048576}')
  NB=$(pgrep -fc 'att_burden_b200\.py')
  NF=$(pgrep -fc 'att_fit_b200\.py')
  NT=$(pgrep -fc 'att_fit_trait\.py')
  NW=$(pgrep -fc 'wd_b200\.sh|wd_fit\.sh|wd_trait\.sh')
  STOPPED=$(ps -u user -o stat=,args= | grep -E 'att_fit|att_burden_b200' \
            | grep -v grep | grep -c '^T')
  echo "$(date -Is) load=$L availGB=$AV swapMB=$SW dfreeGB=$DF wa=$WA% ourRSS=${RSS}G burden=$NB b200fit=$NF trait=$NT wd=$NW stopped=$STOPPED" >> $LOG
  sleep 300
done
