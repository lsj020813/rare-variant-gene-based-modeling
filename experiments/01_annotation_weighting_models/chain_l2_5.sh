#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -u; cd ${PROJECT_ROOT}/work/run_l2; K=${PROJECT_ROOT}/work/run_kcps2fm
while :; do n=$(grep -l POOL_W $K/logs/poolv2_w*.log 2>/dev/null | wc -l); [ "$n" -ge 8 ] && break; sleep 120; done
echo "[$(date '+%m-%d %H:%M')] pool complete: done=$(ls $K/fm/*/*.done | wc -l) skip=$(ls $K/fm/*/*.skip 2>/dev/null | wc -l) fail=$(grep -h 'GATE FAIL' $K/logs/poolv2_w*.log | wc -l)"
( ulimit -v 40000000; python3 l2_5.py > logs/l2_5_full.log 2>&1 ); echo "EXIT=$?" >> logs/l2_5_full.log
grep -q L2_5_DONE logs/l2_5_full.log && touch L2_5_DONE; echo "[$(date '+%m-%d %H:%M')] CHAIN_END"
