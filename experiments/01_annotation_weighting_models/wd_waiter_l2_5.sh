#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
while ! pgrep -f 'python3 l2_5[.]py' >/dev/null; do sleep 20; done
echo "[$(date '+%m-%d %H:%M')] l2_5 detected, starting watchdog" >> ${PROJECT_ROOT}/work/run_l2/logs/wd_waiter.log
exec bash ${PROJECT_ROOT}/work/run_band15/wd5.sh 'l2_5[.]py' ${PROJECT_ROOT}/work/run_l2/wd5 38 40
