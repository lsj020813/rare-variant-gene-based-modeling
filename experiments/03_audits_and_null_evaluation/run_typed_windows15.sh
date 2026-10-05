#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
W=${PROJECT_ROOT}/work/run_meth; PY=${PYTHON:-python3}
echo "START $(date -Is) pid=$$"
nice -n 19 ionice -c3 $PY $W/typed_windows15.py $W/out/typed_windows15.json
rc=$?; echo "END $(date -Is) rc=$rc"; [ $rc -eq 0 ] && echo "typed_windows15 rc=0 $(date -Is)" > $W/out/typed_windows15.done
