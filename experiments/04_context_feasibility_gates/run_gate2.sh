#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
W=${PROJECT_ROOT}/work/gate2; G1=${PROJECT_ROOT}/work/gate1/out
PY=python3
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p $W/out $W/logs
echo "START $(date -Is) pid=$$"
g2_one(){ N=$1; W=${PROJECT_ROOT}/work/gate2; G1=${PROJECT_ROOT}/work/gate1/out
  PY=python3
  [ -s $G1/genes_chr$N.txt ] || return 0
  [ -s $W/out/g2_chr$N.json ] && [ "$N" != "22" ] && { echo "HAVE chr$N"; return 0; }
  nice -n 19 ionice -c3 $PY $W/gate2_context_v2.py $N $G1/genes_chr$N.txt > $W/logs/g2_chr$N.log 2> $W/logs/g2_chr$N.err
  echo "g2 chr$N rc=$? $(tail -1 $W/logs/g2_chr$N.log)"; }
export -f g2_one
seq 22 -1 1 | xargs -P 12 -I{} bash -c 'g2_one {}'
n=$(ls $W/out/g2_chr*.json 2>/dev/null | wc -l)
echo "END $(date -Is) files=$n"
echo "gate2 complete $(date -Is) files=$n" > $W/out/gate2.done
