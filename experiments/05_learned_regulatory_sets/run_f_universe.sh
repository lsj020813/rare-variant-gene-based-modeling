#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
W=${PROJECT_ROOT}/work/fset; mkdir -p $W/out $W/logs
PY=python3
echo "START $(date -Is) pid=$$"
one(){ N=$1; W=${PROJECT_ROOT}/work/fset
  PY=python3
  [ -s $W/out/uni_chr$N.tsv.gz ] && { echo "HAVE chr$N"; return 0; }
  nice -n 19 ionice -c3 $PY $W/f_universe.py $N > $W/logs/uni_chr$N.json 2> $W/logs/uni_chr$N.err
  echo "chr$N rc=$? $(head -c 220 $W/logs/uni_chr$N.json)"; }
export -f one
seq 22 -1 1 | xargs -P 8 -I{} bash -c 'one {}'
n=$(ls $W/out/uni_chr*.tsv.gz 2>/dev/null | wc -l)
echo "END $(date -Is) files=$n"
echo "universe substrate complete $(date -Is) files=$n" > $W/out/universe.done
