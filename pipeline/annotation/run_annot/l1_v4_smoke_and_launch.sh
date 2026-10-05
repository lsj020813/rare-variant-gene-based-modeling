#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -u; cd ${PROJECT_ROOT}/work/run_annot
PY=python3; export OMP_NUM_THREADS=8
$PY -c "import ast;ast.parse(open('l1_train.py').read())" || { echo "GATE FAIL: syntax"; exit 1; }
grep -q 'outer_max=8' l1_train.py || { echo "GATE FAIL: not v4 (no alternating fit)"; exit 1; }
[ "$(pgrep -fc 'l1_trai[n]')" = "0" ] || { echo "GATE FAIL: l1_train still running — stop it first"; exit 1; }
echo "== gradcheck 3 arms (chr22, float64, fixed-β objective) =="
for ARM in spline linear nn; do
  $PY l1_train.py --smoke 22 --arm $ARM --gradcheck-only 2>&1 | grep -E 'GRADCHECK|GATE|Error|Traceback' | sed "s/^/[$ARM] /"
done | tee gc_v4.log
[ "$(grep -c GRADCHECK_ONLY_DONE gc_v4.log)" = "3" ] || { echo "GATE FAIL: gradcheck"; exit 2; }
echo "== smoke spline chr22 (alternating) =="
$PY l1_train.py --smoke 22 --arm spline --maxiter 40 > smoke_v4.log 2>&1
grep -E 'outer|RESULT|Error|GATE|Traceback|L1_DONE' smoke_v4.log | tail -14 | cut -c1-200
grep -q L1_DONE smoke_v4.log || { echo "GATE FAIL: smoke did not finish"; exit 3; }
grep -q "outer 1" smoke_v4.log || { echo "GATE FAIL: optimizer never left outer 0 (still stuck)"; exit 4; }
echo "== launch chain v3 =="
(setsid nohup bash l1_chain.sh >> l1_chain.log 2>&1 </dev/null &); sleep 3
echo "chain=$(pgrep -fc 'bash l1_chai[n].sh')  ($(date '+%F %T'))"
echo LAUNCH_DONE
