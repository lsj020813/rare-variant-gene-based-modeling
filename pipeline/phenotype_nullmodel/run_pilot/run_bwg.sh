#!/usr/bin/env bash
set -uo pipefail
: "${PROJECT_ROOT:?Set PROJECT_ROOT to a generic workspace root}"
PY=${PYTHON:-python3}
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
O="${PROJECT_ROOT}/work/ref/groupfiles_bwg"
build(){ N=$1; [ -s $O/chr$N.B_3kb_re2g.txt ] && { echo "skip chr$N"; return; }
  "$PY" "$SCRIPT_DIR/bwg.py" "$N" && touch $O/chr$N.done || echo "FAIL chr$N"; }
export -f build; export O PY SCRIPT_DIR PROJECT_ROOT
seq 22 -1 1 | xargs -P 6 -I{} bash -c 'build {}'
echo BWG_ALL_DONE
