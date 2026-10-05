#!/bin/bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -u
PAT="$1"; TAG="$2"
ATT=${PROJECT_ROOT}/work/fset/att
LOG=$ATT/logs/wd_$TAG.log
mkdir -p $ATT/logs
ST=0; CRIT=0; IDLE=0
echo "$(date -Is) WD START tag=$TAG pat='$PAT'" >> $LOG
while true; do
  [ -f $ATT/logs/wd_$TAG.stop ] && { echo "$(date -Is) stop marker -> exit" >> $LOG; exit 0; }
  PIDS=$(pgrep -f "$PAT" | tr '\n' ' ')
  if [ -z "$PIDS" ]; then
    IDLE=$((IDLE+1))
    [ "$IDLE" -ge 40 ] && { echo "$(date -Is) 대상 장기 부재 -> exit" >> $LOG; exit 0; }
    sleep 30; continue
  fi
  IDLE=0
  AV=$(free -g | awk '/^Mem:/{print $7}')
  L1=$(cut -d' ' -f1 /proc/loadavg)
  RSS=$(ps -o rss= -p $(echo $PIDS | tr ' ' ',') 2>/dev/null \
        | awk '{s+=$1} END{printf "%d", s/1048576}')
  HI=$(awk -v a="$L1" 'BEGIN{print (a>60)?1:0}')
  LO=$(awk -v a="$AV" 'BEGIN{print (a<60)?1:0}')
  BIG=$(awk -v a="$RSS" 'BEGIN{print (a>60)?1:0}')
  DANGER=$(awk -v a="$AV" 'BEGIN{print (a<25)?1:0}')
  if [ "$DANGER" = "1" ]; then
    CRIT=$((CRIT+1)); echo "$(date -Is) DANGER avail=${AV}G streak=$CRIT" >> $LOG
    if [ "$CRIT" -ge 3 ]; then
      echo "$(date -Is) *** KILL TRIGGERED avail=${AV}G rss=${RSS}G pids=$PIDS" >> $LOG
      for p in $PIDS; do kill -TERM $p 2>/dev/null; done
      sleep 10
      for p in $PIDS; do kill -KILL $p 2>/dev/null; done
      echo "$(date -Is) *** KILL DONE survivors=$(pgrep -f "$PAT" | wc -l)" >> $LOG
      exit 1
    fi
  else CRIT=0; fi
  if [ "$HI" = "1" ] || [ "$LO" = "1" ] || [ "$BIG" = "1" ]; then
    if [ "$ST" = "0" ]; then
      for p in $PIDS; do kill -STOP $p 2>/dev/null; done
      ST=1; echo "$(date -Is) STOP avail=${AV}G rss=${RSS}G load=$L1" >> $LOG; fi
  else
    OKAV=$(awk -v a="$AV" 'BEGIN{print (a>80)?1:0}')
    OKL=$(awk -v a="$L1" 'BEGIN{print (a<45)?1:0}')
    if [ "$ST" = "1" ] && [ "$OKAV" = "1" ] && [ "$OKL" = "1" ]; then
      for p in $PIDS; do kill -CONT $p 2>/dev/null; done
      ST=0; echo "$(date -Is) CONT avail=${AV}G rss=${RSS}G load=$L1" >> $LOG; fi
  fi
  echo "$(date -Is) avail=${AV}G rss=${RSS}G load1=$L1 n=$(echo $PIDS|wc -w) stopped=$ST" >> $LOG
  sleep 30
done
