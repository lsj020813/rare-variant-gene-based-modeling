#!/usr/bin/env bash
set -uo pipefail
: "${PROJECT_ROOT:?Set PROJECT_ROOT to a generic workspace root}"
RUN_DIR=${RUN_DIR:-"${PROJECT_ROOT}/work/run_pilot"}
mkdir -p "$RUN_DIR"
exec 9>"$RUN_DIR/.pilot.lock"
flock -n 9 || { echo "another instance running"; exit 0; }

R=${REFERENCE_DIR:-"${PROJECT_ROOT}/work/ref"}
OUT=$R/saige_pilot
P=$R/groupfiles_pilot
LOG="$RUN_DIR"
mkdir -p "$OUT"
SAIGE_RUNNER=${SAIGE_RUNNER:-Rscript}
SAIGE_STEP2_SCRIPT=${SAIGE_STEP2_SCRIPT:?Set SAIGE_STEP2_SCRIPT}

run_one() {
  local TAG=$1 TRAIT=$2 CH=$3 GF=$4
  local O="$OUT/$TAG"
  [ -s "$O.done" ] && { echo "skip $TAG"; return 0; }
  "$SAIGE_RUNNER" "$SAIGE_STEP2_SCRIPT" \
    --vcfFile="$R/band_vcf/chr$CH.band.vcf.gz" \
    --vcfFileIndex="$R/band_vcf/chr$CH.band.vcf.gz.csi" \
    --vcfField=DS \
    --chrom="$CH" \
    --GMMATmodelFile="$R/saige_step1/${TRAIT}_v2.rda" \
    --varianceRatioFile="$R/saige_step1/${TRAIT}_v2.varianceRatio.txt" \
    --sparseGRMFile="${GRM_FILE:?Set GRM_FILE}" \
    --sparseGRMSampleIDFile="${GRM_SAMPLE_IDS_FILE:?Set GRM_SAMPLE_IDS_FILE}" \
    --groupFile="$GF" \
    --annotation_in_groupTest=noncoding \
    --maxMAF_in_groupTest=0.01 \
    --is_output_moreDetails=TRUE \
    --LOCO=FALSE \
    --SAIGEOutputFile="$O" > "$O.log" 2>&1
  local rc=$?
  local n=$(grep -vc '^Region' "$O" 2>/dev/null || echo 0)
  if [ "$rc" -ne 0 ] || [ "${n:-0}" -lt 1 ]; then
    echo "FAIL $TAG rc=$rc rows=$n"; return 1
  fi
  : > "$O.done"; echo "OK $TAG rows=$n"
}
export -f run_one; export R OUT SAIGE_RUNNER SAIGE_STEP2_SCRIPT GRM_FILE GRM_SAMPLE_IDS_FILE

: > "$LOG/worklist.txt"
for ARM in B_3kb_re2g C_nearest D_3kb_only; do
  echo "pilot_${ARM}_tchl_chr19|tchl|19|$P/chr19.${ARM}.txt" >> "$LOG/worklist.txt"
done
for PART in 001 002 003 004 005; do
  echo "recover_tchl_chr2_part${PART}|tchl|2|$R/groupfiles_chunks/chr2.part${PART}.txt" >> "$LOG/worklist.txt"
done
wc -l < "$LOG/worklist.txt"

POOL=${POOL:-8}
tr '|' ' ' < "$LOG/worklist.txt" | xargs -P "$POOL" -L1 bash -c 'run_one "$@"' _
echo "PILOT_RUN_COMPLETE done=$(ls $OUT/*.done 2>/dev/null | wc -l)/8"
