#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -u
cd ${PROJECT_ROOT}/work/run_annot_bbj
E6=${PROJECT_ROOT}/work/run_band15/model_v10_out/e6_runs; L1=${PROJECT_ROOT}/work/ref/annot/bbj_l1
LOG=chain_full.log; log(){ echo "$(date '+%F %T') $*" >> $LOG; }
log "CHAIN8 START (wait LC -> D3 subset -> sens_phi17 -> sens_ecdfeval)"
while [ ! -s $E6/primary_f0_lc/rc.txt ]; do sleep 30; done
log "primary_lc rc=$(cat $E6/primary_f0_lc/rc.txt)"; python3 summarize_results.py > /dev/null 2>&1
FO='{"0":["13"],"1":["16"],"2":["7"],"3":["19"],"4":[]}'
bash v10_run.sh d3_sub_f0 _d3 --fold 0 --inner-va-folds 1 --folds-override "$FO" --phi-columns v9-eighteen --c1-ecdf-reference training --train-bands all \
    --simulate --oracle-spec $E6/d3_spec_primary_f0.json --effect-grid 0 0.25 0.5 1 2 --simulation-reps 20
log "d3_sub rc=$(cat $E6/d3_sub_f0/rc.txt)"
bash v10_run.sh sens_phi17_f0 "" --fold 0 --phi-columns v8-seventeen --c1-ecdf-reference training --train-bands all; log "sens_phi17 rc=$(cat $E6/sens_phi17_f0/rc.txt)"
python3 summarize_results.py > /dev/null 2>&1
bash v10_run.sh sens_ecdfeval_f0 "" --fold 0 --phi-columns v9-eighteen --c1-ecdf-reference evaluation --train-bands all; log "sens_ecdfeval rc=$(cat $E6/sens_ecdfeval_f0/rc.txt)"
python3 summarize_results.py > /dev/null 2>&1
log "CHAIN8 END"
