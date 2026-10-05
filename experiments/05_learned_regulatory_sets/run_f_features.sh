#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
S=$1
W=${PROJECT_ROOT}/work/fset; mkdir -p $W/out $W/logs
PY=python3
bash $W/f_remap_split.sh
echo "START $(date -Is) pid=$$ sample=$S"
one(){ N=$1; S=$2; W=${PROJECT_ROOT}/work/fset
  PY=python3
  [ -s $W/out/feat_chr$N.tsv.gz ] && { echo "HAVE chr$N"; return 0; }
  nice -n 19 ionice -c3 $PY $W/f_features_v2.py $N $S > $W/logs/feat_chr$N.json 2> $W/logs/feat_chr$N.err
  echo "feat chr$N rc=$? $(head -c 260 $W/logs/feat_chr$N.json)"; }
export -f one
cut -f2 $S | tail -n+2 | sort -u | xargs -P 8 -I{} bash -c "one {} $S"
n=$(ls $W/out/feat_chr*.tsv.gz 2>/dev/null | wc -l)
echo "END $(date -Is) files=$n"
echo "features complete $(date -Is) files=$n" > $W/out/features.done
