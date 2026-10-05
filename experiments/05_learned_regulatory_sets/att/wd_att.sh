#!/bin/bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
ATT=${PROJECT_ROOT}/work/fset/att
LOG=$ATT/logs/watchdog.log
PAT='[a]tt_burden.py'
ST=0
while true; do
  if [ -f $ATT/burden_stage.done ]; then
    echo "$(date -Is) done marker -> watchdog exit" >> $LOG; exit 0; fi
  L1=$(cut -d' ' -f1 /proc/loadavg)
  AV=$(free -g | awk '/^Mem:/{print $7}')
  PIDS=$(pgrep -f "$PAT" | tr '\n' ' ')
  HI=$(awk -v a="$L1" 'BEGIN{print (a>40)?1:0}')
  LO=$(awk -v a="$AV" 'BEGIN{print (a<20)?1:0}')
  if [ -n "$PIDS" ]; then
    if [ "$HI" = "1" ] || [ "$LO" = "1" ]; then
      if [ "$ST" = "0" ]; then kill -STOP $PIDS 2>/dev/null; ST=1;
        echo "$(date -Is) STOP load=$L1 avail=${AV}G" >> $LOG; fi
    else
      if [ "$ST" = "1" ]; then kill -CONT $PIDS 2>/dev/null; ST=0;
        echo "$(date -Is) CONT load=$L1 avail=${AV}G" >> $LOG; fi
    fi
  fi
  echo "$(date -Is) load1=$L1 availGB=$AV n=$(echo $PIDS|wc -w) stopped=$ST" >> $LOG
  sleep 60
done
