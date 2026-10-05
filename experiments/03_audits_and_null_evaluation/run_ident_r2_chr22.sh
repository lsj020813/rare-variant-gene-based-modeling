#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
W=${PROJECT_ROOT}/work/run_meth; mkdir -p $W/out
echo "START $(date -Is) host=$(hostname) pid=$$"
nice -n 10 python $W/ident_r2_chr22.py $W/out/ident_r2_chr22.json
rc=$?; echo "END $(date -Is) rc=$rc"
[ $rc -eq 0 ] && echo "ident_r2_chr22 rc=0 $(date -Is)" > $W/out/ident_r2_chr22.done
