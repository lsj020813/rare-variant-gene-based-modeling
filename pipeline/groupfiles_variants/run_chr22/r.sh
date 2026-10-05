#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
exec 9>${PROJECT_ROOT}/work/run_chr22/.lock
flock -n 9 || { echo "another instance running"; exit 0; }

R=${PROJECT_ROOT}/work/ref
S1=$R/saige_step1_v4
CH=$R/groupfiles_chunks
BAND=$R/band_vcf
OUT=$R/saige_step2_v4
mkdir -p "$OUT"
export UDOCKER_DIR=${CONTAINER_STATE_DIR:?Set CONTAINER_STATE_DIR}
UD="PROOT_NO_SECCOMP=1 ${UDOCKER_BIN:-udocker} run --volume=/data:/data ${SAIGE_IMAGE:?Set SAIGE_IMAGE}"
POOL=${POOL:-8}

for T in htn dm lip tchl; do
  for f in "$S1/${T}_v4.rda" "$S1/${T}_v4.varianceRatio.txt"; do
    [ -s "$f" ] || { echo "PREFLIGHT FAIL missing $f"; exit 3; }
  done
done
[ -s "$BAND/chr22.band.vcf.gz.csi" ] || { echo "PREFLIGHT FAIL missing chr22 index"; exit 3; }
echo "preflight OK"

run_one() {
  local T=$1 P=$2
  local G=$CH/chr22.$P.txt
  local O=$OUT/$T.chr22.$P
  [ -s "$O.done" ] && { echo "  skip $T.$P"; return 0; }
  eval $UD step2_SPAtests.R \
    --vcfFile="$BAND/chr22.band.vcf.gz" --vcfFileIndex="$BAND/chr22.band.vcf.gz.csi" \
    --vcfField="DS" --chrom="22" --AlleleOrder=ref-first \
    --minMAF=0 --minMAC=0.5 --LOCO=FALSE --is_fastTest=TRUE \
    --GMMATmodelFile="$S1/${T}_v4.rda" --varianceRatioFile="$S1/${T}_v4.varianceRatio.txt" \
    --groupFile="$G" --annotation_in_groupTest="all" --maxMAF_in_groupTest=0.01 \
    --is_output_markerList_in_groupTest=TRUE \
    --SAIGEOutputFile="$O" > "$O.log" 2>&1
  local rc=$?
  local nrow=0; [ -s "$O" ] && nrow=$(( $(wc -l < "$O") - 1 ))
  if [ "$rc" -ne 0 ] || [ "$nrow" -le 0 ]; then
    echo "  FAIL $T.$P rc=$rc rows=$nrow"; tail -4 "$O.log"; return 1
  fi
  touch "$O.done"; echo "  OK $T.$P rows=$nrow"
}
export -f run_one; export CH BAND OUT S1 UD
for T in htn dm lip tchl; do for P in part001 part002; do echo "$T|$P"; done; done \
  | tr '|' ' ' | xargs -P "$POOL" -L1 bash -c 'run_one "$@"' _

n=$(ls $OUT/*.done 2>/dev/null | wc -l)
echo "completed: $n/8"
[ "$n" -eq 8 ] || { echo "GATE FAIL: $n/8"; exit 6; }
echo "CHR22_V4_COMPLETE"
