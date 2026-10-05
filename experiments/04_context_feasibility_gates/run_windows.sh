#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
W=${PROJECT_ROOT}/work/gate1; mkdir -p $W/out $W/logs
PY=python3
echo "START $(date -Is) pid=$$"
run_one() {
  N=$1; W=${PROJECT_ROOT}/work/gate1
  PY=python3
  [ -s "$W/out/window_summary_chr$N.tsv" ] && { echo "HAVE chr$N"; return 0; }
  nice -n 19 ionice -c3 $PY $W/gate1_windows.py $N > $W/logs/windows_chr$N.json 2> $W/logs/windows_chr$N.err
  rc=$?
  echo "chr$N rc=$rc $(cat $W/logs/windows_chr$N.json 2>/dev/null | head -c 300)"
  [ $rc -ne 0 ] && head -3 $W/logs/windows_chr$N.err
  return 0
}
export -f run_one
seq 22 -1 1 | xargs -P 4 -I{} bash -c 'run_one {}'
echo "END $(date -Is)"
n=$(ls $W/out/window_summary_chr*.tsv 2>/dev/null | wc -l)
echo "windows built for $n chromosomes $(date -Is)" > $W/out/windows.done
