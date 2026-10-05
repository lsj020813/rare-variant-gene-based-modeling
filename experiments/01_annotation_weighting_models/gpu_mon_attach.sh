#!/bin/bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
n=0
while ! pgrep -f 'enformer_score_v2[.]py' >/dev/null; do sleep 2; n=$((n+1)); if [ $n -ge 90 ]; then echo "WARN gpu_mon target never appeared" >> ${PROJECT_ROOT}/work/run_trackB/gpu_mon.log; exit 2; fi; done
exec bash ${PROJECT_ROOT}/work/run_trackB/gpu_mon.sh
