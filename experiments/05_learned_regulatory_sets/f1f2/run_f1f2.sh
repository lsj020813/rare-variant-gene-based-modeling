#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -o pipefail
export TMPDIR=${PROJECT_ROOT}/work/tmp TMP=${PROJECT_ROOT}/work/tmp TEMP=${PROJECT_ROOT}/work/tmp
W=${PROJECT_ROOT}/work/fset/f1f2; L=$W/logs; mkdir -p $L
PY=python3
export F_WORKERS=4 F_NQ=10000
ulimit -v $((24*1024*1024))
run() { local tag=$1; shift; echo "$(date -Is) START $tag" >> $L/run.log; nice -n 19 ionice -c3 $PY "$@" > $L/$tag.log 2>&1; echo "$(date -Is) END $tag exit=$?" >> $L/run.log; }
F_PERDOM=100 F_BLOCKNULL=1 run f1U_ub_inub $W/f1_tendency.py ub_inub U
run f1C_ub_inub $W/f1_tendency.py ub_inub C
run f2U_ub_inub $W/f2_modelU.py ub_inub
run f2C_ub_inub $W/f2_modelC.py ub_inub
F_PERDOM=50 F_BLOCKNULL=0 run f1U_ub      $W/f1_tendency.py ub U
run f1C_ub      $W/f1_tendency.py ub C
F_PERDOM=50 F_BLOCKNULL=0 run f1U_ua      $W/f1_tendency.py ua U
run f1C_ua      $W/f1_tendency.py ua C
run f2U_ub      $W/f2_modelU.py ub sens
run f2C_ub      $W/f2_modelC.py ub sens
run f2U_ua      $W/f2_modelU.py ua sens
run f2C_ua      $W/f2_modelC.py ua sens
touch $W/ALL.done; echo "$(date -Is) ALL DONE" >> $L/run.log
