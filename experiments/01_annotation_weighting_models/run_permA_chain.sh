#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -o pipefail
B=${PROJECT_ROOT}/work/run_l3b
for n in A1 A2 A3; do
  echo "=== perm$n dispatch $(date '+%m-%d %H:%M:%S') load=$(cut -d' ' -f1 /proc/loadavg) ==="
  rm -f "$B/out/STOP"
  bash "$B/run_arm_seq2.sh" "perm$n" 16 45 70
  rc=$?
  echo "=== perm$n rc=$rc $(date '+%m-%d %H:%M:%S') ==="
  [ $rc -ne 0 ] && { echo "CHAIN ABORTED at perm$n"; exit $rc; }
done
echo "PERMA_CHAIN_DONE $(date '+%m-%d %H:%M:%S')"
