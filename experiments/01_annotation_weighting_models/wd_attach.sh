#!/bin/bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
PAT="$1"; WD="$2"; LOAD="$3"; RSS="$4"; n=0
while ! pgrep -f "$PAT" >/dev/null; do sleep 2; n=$((n+1)); if [ $n -ge 90 ]; then echo "WARN watchdog target never appeared: $PAT" >> "$WD/watchdog.log"; exit 2; fi; done
exec bash ${PROJECT_ROOT}/work/run_band15/wd5.sh "$PAT" "$WD" "$LOAD" "$RSS"
