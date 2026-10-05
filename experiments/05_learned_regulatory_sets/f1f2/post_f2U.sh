#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
export TMPDIR=${PROJECT_ROOT}/work/tmp TMP=${PROJECT_ROOT}/work/tmp TEMP=${PROJECT_ROOT}/work/tmp
W=${PROJECT_ROOT}/work/fset/f1f2; L=$W/logs; PY=python3; export F_WORKERS=4
while [ ! -f $W/ALL.done ]; do sleep 60; done
mkdir -p $W/ub_inub/f2U_attempt1_connected_components && cp $W/ub_inub/f2U_consensus_sizes.csv $W/ub_inub/f2U_consensus_assign.npz $W/ub_inub/f2U_attempt1_connected_components/ 2>/dev/null
echo "$(date -Is) START f2U_ub_inub_rerun" >> $L/run.log
nice -n 19 ionice -c3 $PY $W/f2_modelU.py ub_inub > $L/f2U_ub_inub_rerun.log 2>&1; echo "$(date -Is) END f2U_ub_inub_rerun exit=$?" >> $L/run.log
touch $W/POST.done
