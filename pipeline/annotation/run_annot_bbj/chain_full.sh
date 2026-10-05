#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -u
cd ${PROJECT_ROOT}/work/run_annot_bbj
PY=python3
A=${PROJECT_ROOT}/work/ref/annot; L1=$A/bbj_l1; E6=${PROJECT_ROOT}/work/run_band15/model_v10_out/e6_runs
LOG=chain_full.log; log(){ echo "$(date '+%F %T') $*" >> $LOG; }
log "CHAIN START"
for i in $(seq 1 540); do
  a=$(ls $A/extract_bbj/chr{1..20}.annot.tsv.done 2>/dev/null | wc -l); b=$(ls $A/extract_bbj/chr{1..20}.cadd.tsv.done 2>/dev/null | wc -l)
  g=$(ls $A/gpn_bbj/chr{1..20}.gpn.tsv.done 2>/dev/null | wc -l); t=$(ls $A/t2_bbj/chr{1..20}.t2.tsv.done 2>/dev/null | wc -l)
  [ "$a" -eq 20 ] && [ "$b" -eq 20 ] && [ "$g" -eq 20 ] && [ "$t" -eq 20 ] && break; sleep 20; done
log "inputs: favor $a/$b gpn $g t2 $t"
[ "$a" -eq 20 ] && [ "$g" -eq 20 ] && [ "$t" -eq 20 ] || { log "GATE FAIL: inputs incomplete"; exit 2; }
( ulimit -v 40000000; $PY assemble_bbj.py --chrs 20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1 --out bbj_annot.tsv.gz > asm_full.log 2>&1 ); log "assemble rc=$?"
[ -s $L1/bbj_annot.tsv.gz.done ] || { log "GATE FAIL: assemble"; exit 3; }
( ulimit -v 16000000; $PY pip_bbj_B.py --annot bbj_annot.tsv.gz --outdir bbj_pip --traits-out traits.tsv > pipB_full.log 2>&1 ); log "pipB rc=$?"
[ -s $L1/bbj_pip/.done ] || { log "GATE FAIL: pipB"; exit 4; }
echo "$(date '+%F %T') INPUTS_READY" >> $LOG
bash v10_run.sh primary_f0 "" --fold 0 --phi-columns v9-eighteen --c1-ecdf-reference training --train-bands all --learning-curves
log "primary rc=$(cat $E6/primary_f0/rc.txt 2>/dev/null)"
[ -s $E6/primary_f0/L1_DONE ] || { log "GATE FAIL: primary did not finish — D3 skipped; sensitivity continues"; }
if [ -s $E6/primary_f0/L1_DONE ]; then
  $PY d3_spec.py $E6/primary_f0/phi.json $L1/bbj_pip $L1/traits.tsv $E6/d3_spec_primary_f0.json > d3_spec.log 2>&1 && log "d3_spec ok: $(cat d3_spec.log)" || log "d3_spec FAIL"
  if [ -s $E6/d3_spec_primary_f0.json ]; then
    ( bash v10_run.sh d3_f0 "" --fold 0 --phi-columns v9-eighteen --c1-ecdf-reference training --train-bands all \
        --simulate --oracle-spec $E6/d3_spec_primary_f0.json --effect-grid 0 0.25 0.5 1 2 --simulation-reps 50; log "d3 rc=$(cat $E6/d3_f0/rc.txt)" ) &
  fi
fi
bash v10_run.sh sens_phi17_f0 "" --fold 0 --phi-columns v8-seventeen --c1-ecdf-reference training --train-bands all;   log "sens_phi17 rc=$(cat $E6/sens_phi17_f0/rc.txt)"
bash v10_run.sh sens_ecdfeval_f0 "" --fold 0 --phi-columns v9-eighteen --c1-ecdf-reference evaluation --train-bands all; log "sens_ecdfeval rc=$(cat $E6/sens_ecdfeval_f0/rc.txt)"
bash v10_run.sh sens_tb15_f0 "" --fold 0 --phi-columns v9-eighteen --c1-ecdf-reference training --train-bands 1-5;     log "sens_tb15 rc=$(cat $E6/sens_tb15_f0/rc.txt)"
bash v10_run.sh sens_tb011_f0 "" --fold 0 --phi-columns v9-eighteen --c1-ecdf-reference training --train-bands 0.1-1;  log "sens_tb011 rc=$(cat $E6/sens_tb011_f0/rc.txt)"
wait
log "CHAIN END"
