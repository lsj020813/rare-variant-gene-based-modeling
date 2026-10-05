#!/bin/bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
ATT=${PROJECT_ROOT}/work/fset/att
OUT=$ATT/b200
LOG=$OUT/logs/watchdog.log
PAT='[a]tt_burden_b200.py'
ST=0
while true; do
  if [ -f $OUT/burden_stage.done ]; then
    echo "$(date -Is) done marker -> watchdog exit" >> $LOG; exit 0; fi
  L1=$(cut -d' ' -f1 /proc/loadavg)
  AV=$(free -g | awk '/^Mem:/{print $7}')
  DF=$(df -BG --output=avail /data | tail -1 | tr -dc '0-9')
  PIDS=$(pgrep -f "$PAT" | tr '\n' ' ')
  HI=$(awk -v a="$L1" 'BEGIN{print (a>40)?1:0}')
  LO=$(awk -v a="$AV" 'BEGIN{print (a<20)?1:0}')
  DL=$(awk -v a="$DF" 'BEGIN{print (a<500)?1:0}')
  if [ -n "$PIDS" ]; then
    if [ "$HI" = "1" ] || [ "$LO" = "1" ] || [ "$DL" = "1" ]; then
      if [ "$ST" = "0" ]; then kill -STOP $PIDS 2>/dev/null; ST=1;
        echo "$(date -Is) STOP load=$L1 avail=${AV}G dfree=${DF}G" >> $LOG; fi
    else
      if [ "$ST" = "1" ]; then kill -CONT $PIDS 2>/dev/null; ST=0;
        echo "$(date -Is) CONT load=$L1 avail=${AV}G dfree=${DF}G" >> $LOG; fi
    fi
  fi
  echo "$(date -Is) load1=$L1 availGB=$AV dfreeGB=$DF n=$(echo $PIDS|wc -w) stopped=$ST" >> $LOG
  sleep 60
done
