#!/bin/bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
while pgrep -f 'enformer_score_v2[.]py' >/dev/null; do echo "$(date +%T) $(nvidia-smi --query-gpu=memory.used,utilization.gpu,temperature.gpu --format=csv,noheader)"; sleep 5; done >> ${PROJECT_ROOT}/work/run_trackB/gpu_mon.log
