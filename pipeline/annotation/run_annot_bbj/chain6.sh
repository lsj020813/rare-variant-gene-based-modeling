#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -u
cd ${PROJECT_ROOT}/work/run_annot_bbj
PY=python3
A=${PROJECT_ROOT}/work/ref/annot; L1=$A/bbj_l1; E6=${PROJECT_ROOT}/work/run_band15/model_v10_out/e6_runs
LOG=chain_full.log; log(){ echo "$(date '+%F %T') $*" >> $LOG; }
log "CHAIN6 START (primary rerun after E6C memory patch)"
mv $E6/primary_f0 $E6/primary_f0_killed_0445 2>/dev/null; rm -rf $E6/sens_phi17_f0
bash v10_run.sh primary_f0 "" --fold 0 --phi-columns v9-eighteen --c1-ecdf-reference training --train-bands all
log "primary rc=$(cat $E6/primary_f0/rc.txt 2>/dev/null)"
if [ -s $E6/primary_f0/L1_DONE ]; then
  $PY d3_spec.py $E6/primary_f0/phi.json $L1/bbj_pip $L1/traits.tsv $E6/d3_spec_primary_f0.json > d3_spec.log 2>&1 && log "d3_spec ok: $(cat d3_spec.log)" || log "d3_spec FAIL"
  if [ -s $E6/d3_spec_primary_f0.json ]; then
    ( bash v10_run.sh d3_f0 "" --fold 0 --phi-columns v9-eighteen --c1-ecdf-reference training --train-bands all \
        --simulate --oracle-spec $E6/d3_spec_primary_f0.json --effect-grid 0 0.25 0.5 1 2 --simulation-reps 50; log "d3 rc=$(cat $E6/d3_f0/rc.txt)" ) &
  fi
else
  log "GATE FAIL: primary did not finish — D3 skipped; sensitivity continues"
fi
bash v10_run.sh primary_f0_lc "" --fold 0 --phi-columns v9-eighteen --c1-ecdf-reference training --train-bands all --learning-curves; log "primary_lc rc=$(cat $E6/primary_f0_lc/rc.txt)"
bash v10_run.sh sens_phi17_f0 "" --fold 0 --phi-columns v8-seventeen --c1-ecdf-reference training --train-bands all;   log "sens_phi17 rc=$(cat $E6/sens_phi17_f0/rc.txt)"
bash v10_run.sh sens_ecdfeval_f0 "" --fold 0 --phi-columns v9-eighteen --c1-ecdf-reference evaluation --train-bands all; log "sens_ecdfeval rc=$(cat $E6/sens_ecdfeval_f0/rc.txt)"
bash v10_run.sh sens_tb15_f0 "" --fold 0 --phi-columns v9-eighteen --c1-ecdf-reference training --train-bands 1-5;     log "sens_tb15 rc=$(cat $E6/sens_tb15_f0/rc.txt)"
bash v10_run.sh sens_tb011_f0 "" --fold 0 --phi-columns v9-eighteen --c1-ecdf-reference training --train-bands 0.1-1;  log "sens_tb011 rc=$(cat $E6/sens_tb011_f0/rc.txt)"
wait
log "CHAIN6 END"
