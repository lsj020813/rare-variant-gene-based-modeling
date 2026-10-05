#!/usr/bin/env bash
set -uo pipefail
: "${PROJECT_ROOT:?Set PROJECT_ROOT to a generic workspace root}"
RUN_DIR=${RUN_DIR:-"${PROJECT_ROOT}/work/run_pilot"}
mkdir -p "$RUN_DIR"
exec 9>"$RUN_DIR/.lock2"
flock -n 9 || { echo "another instance running"; exit 0; }

R=${REFERENCE_DIR:-"${PROJECT_ROOT}/work/ref"}
S1=$R/saige_step1_v4
GP=$R/groupfiles_pilot
BAND=$R/band_vcf
OUT=$R/saige_step2_pilot
mkdir -p "$OUT"
SAIGE_RUNNER=${SAIGE_RUNNER:-Rscript}
SAIGE_STEP2_SCRIPT=${SAIGE_STEP2_SCRIPT:?Set SAIGE_STEP2_SCRIPT}
POOL=${POOL:-6}
ARMS=${ARMS:-"B_3kb_re2g C_nearest D_3kb_only"}
T=tchl

for f in "$S1/${T}_v4.rda" "$S1/${T}_v4.varianceRatio.txt" \
         "$BAND/chr19.band.vcf.gz" "$BAND/chr19.band.vcf.gz.csi"; do
  [ -s "$f" ] || { echo "PREFLIGHT FAIL missing $f"; exit 3; }
done
for A in $ARMS; do
  [ -s "$GP/chr19.$A.txt" ] || { echo "PREFLIGHT FAIL missing $GP/chr19.$A.txt"; exit 3; }
done
echo "preflight OK  arms=$ARMS  pool=$POOL"

run_arm() {
  local A=$1
  local G=$GP/chr19.$A.txt
  local O=$OUT/tchl.chr19.$A
  [ -s "$O.done" ] && { echo "  skip $A"; return 0; }
  "$SAIGE_RUNNER" "$SAIGE_STEP2_SCRIPT" \
    --vcfFile="$BAND/chr19.band.vcf.gz" --vcfFileIndex="$BAND/chr19.band.vcf.gz.csi" \
    --vcfField="DS" --chrom="19" --AlleleOrder=ref-first \
    --minMAF=0 --minMAC=0.5 --LOCO=FALSE --is_fastTest=TRUE \
    --GMMATmodelFile="$S1/tchl_v4.rda" --varianceRatioFile="$S1/tchl_v4.varianceRatio.txt" \
    --groupFile="$G" --annotation_in_groupTest="all" --maxMAF_in_groupTest=0.01 \
    --is_output_markerList_in_groupTest=TRUE \
    --SAIGEOutputFile="$O" > "$O.log" 2>&1
  local rc=$?
  local exp=$(awk '$2=="var"{print $1}' "$G" | sort -u | wc -l)
  local nrow=0; [ -s "$O" ] && nrow=$(( $(wc -l < "$O") - 1 ))
  if [ "$rc" -ne 0 ] || [ "$nrow" -le 0 ]; then
    echo "  FAIL $A rc=$rc rows=$nrow/$exp"; tail -3 "$O.log"; return 1
  fi
  touch "$O.done"
  [ "$nrow" -eq "$exp" ] && echo "  OK $A rows=$nrow" || echo "  PARTIAL $A $nrow/$exp"
}
export -f run_arm; export GP BAND OUT S1 SAIGE_RUNNER SAIGE_STEP2_SCRIPT
echo "$ARMS" | tr ' ' '\n' | xargs -P "$POOL" -I{} bash -c 'run_arm {}'

n=$(ls $OUT/*.done 2>/dev/null | wc -l)
echo "completed: $n"
echo "PILOT_RUN_COMPLETE"
