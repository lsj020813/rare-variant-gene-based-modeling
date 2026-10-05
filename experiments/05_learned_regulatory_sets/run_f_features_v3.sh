#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
S=$1
W=${PROJECT_ROOT}/work/fset; O=$W/primary; mkdir -p $W/logs
PY=python3
[ -s $W/out/remap_by_chr/chr22.bed ] || bash $W/f_remap_split.sh
echo "START $(date -Is) pid=$$ sample=$S load=$(cut -d' ' -f1 /proc/loadavg)"
one(){ N=$1; S=$2; W=${PROJECT_ROOT}/work/fset; O=$W/primary
  PY=python3
  [ -s $O/feat_chr$N.tsv.gz ] && [ -s $O/vset_re2g_chr$N.tsv.gz ] && { echo "HAVE chr$N"; return 0; }
  nice -n 19 ionice -c3 $PY $W/f_features_v3.py $N $S > $W/logs/pfeat_chr$N.json 2> $W/logs/pfeat_chr$N.err
  echo "feat chr$N rc=$? $(date -Is) $(head -c 300 $W/logs/pfeat_chr$N.json)"; }
export -f one
cut -f2 $S | tail -n+2 | sort -u | xargs -P 8 -I{} bash -c "one {} $S"
n=$(ls $O/feat_chr*.tsv.gz 2>/dev/null | wc -l)
echo "END $(date -Is) files=$n"
echo "primary features complete $(date -Is) files=$n" > $O/features.done
