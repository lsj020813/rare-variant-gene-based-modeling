#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -u
cd ${PROJECT_ROOT}/work/run_annot_bbj
E6=${PROJECT_ROOT}/work/run_band15/model_v10_out/e6_runs
LOG=chain_full.log; log(){ echo "$(date '+%F %T') $*" >> $LOG; }
log "CHAIN7 START (ram cap 60G per job; LC -> sens_phi17 -> sens_ecdfeval; no D3)"
for d in primary_f0_lc sens_phi17_f0 sens_ecdfeval_f0; do [ -d $E6/$d ] && [ ! -s $E6/$d/L1_DONE ] && mv $E6/$d $E6/${d}_killed38G; done
bash v10_run.sh primary_f0_lc "" --fold 0 --phi-columns v9-eighteen --c1-ecdf-reference training --train-bands all --learning-curves; log "primary_lc rc=$(cat $E6/primary_f0_lc/rc.txt)"
python3 summarize_results.py > /dev/null 2>&1
bash v10_run.sh sens_phi17_f0 "" --fold 0 --phi-columns v8-seventeen --c1-ecdf-reference training --train-bands all; log "sens_phi17 rc=$(cat $E6/sens_phi17_f0/rc.txt)"
python3 summarize_results.py > /dev/null 2>&1
bash v10_run.sh sens_ecdfeval_f0 "" --fold 0 --phi-columns v9-eighteen --c1-ecdf-reference evaluation --train-bands all; log "sens_ecdfeval rc=$(cat $E6/sens_ecdfeval_f0/rc.txt)"
python3 summarize_results.py > /dev/null 2>&1
log "CHAIN7 END"
