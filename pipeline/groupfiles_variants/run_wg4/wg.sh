#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
exec 9>${PROJECT_ROOT}/work/run_wg4/.lock
flock -n 9 || { echo "another instance running"; exit 0; }

R=${PROJECT_ROOT}/work/ref
S1=$R/saige_step1_v4
CH=$R/groupfiles_chunks
BAND=$R/band_vcf
OUT=$R/saige_step2_v4
mkdir -p "$OUT"
export UDOCKER_DIR=${CONTAINER_STATE_DIR:?Set CONTAINER_STATE_DIR}
UD="PROOT_NO_SECCOMP=1 ${UDOCKER_BIN:-udocker} run --volume=/data:/data ${SAIGE_IMAGE:?Set SAIGE_IMAGE}"
POOL=${POOL:-64}

for T in htn dm lip tchl; do
  for f in "$S1/${T}_v4.rda" "$S1/${T}_v4.varianceRatio.txt"; do
    [ -s "$f" ] || { echo "PREFLIGHT FAIL missing $f"; exit 3; }
  done
done
for N in $(seq 1 22); do
  [ -s "$BAND/chr$N.band.vcf.gz" ] && [ -s "$BAND/chr$N.band.vcf.gz.csi" ] \
    || { echo "PREFLIGHT FAIL band chr$N"; exit 3; }
  ls $CH/chr$N.part*.txt >/dev/null 2>&1 || { echo "PREFLIGHT FAIL groupfile chr$N"; exit 3; }
done
NCH=$(ls $CH/chr*.part*.txt | wc -l)
echo "preflight OK  chunks=$NCH  expected runs=$(( NCH * 4 ))  pool=$POOL"

run_one() {
  local T=$1 N=$2 P=$3
  local G=$CH/chr$N.$P.txt
  local O=$OUT/$T.chr$N.$P
  [ -s "$O.done" ] && return 0
  eval $UD step2_SPAtests.R \
    --vcfFile="$BAND/chr$N.band.vcf.gz" --vcfFileIndex="$BAND/chr$N.band.vcf.gz.csi" \
    --vcfField="DS" --chrom="$N" --AlleleOrder=ref-first \
    --minMAF=0 --minMAC=0.5 --LOCO=FALSE --is_fastTest=TRUE \
    --GMMATmodelFile="$S1/${T}_v4.rda" --varianceRatioFile="$S1/${T}_v4.varianceRatio.txt" \
    --groupFile="$G" --annotation_in_groupTest="all" --maxMAF_in_groupTest=0.01 \
    --is_output_markerList_in_groupTest=TRUE \
    --SAIGEOutputFile="$O" > "$O.log" 2>&1
  local rc=$?
  local exp=$(awk '$2=="var"{print $1}' "$G" | sort -u | wc -l)
  local nrow=0; [ -s "$O" ] && nrow=$(( $(wc -l < "$O") - 1 ))
  if [ "$rc" -ne 0 ] || [ "$nrow" -le 0 ]; then
    echo "  FAIL $T.chr$N.$P rc=$rc rows=$nrow/$exp"; tail -3 "$O.log"; return 1
  fi
  touch "$O.done"
  [ "$nrow" -eq "$exp" ] && echo "  OK $T.chr$N.$P $nrow" \
                         || echo "  PARTIAL $T.chr$N.$P $nrow/$exp"
}
export -f run_one; export CH BAND OUT S1 UD

: > ${PROJECT_ROOT}/work/run_wg4/worklist.txt
for N in $(seq 1 22); do
  for f in $CH/chr$N.part*.txt; do
    P=$(basename "$f" .txt); P=${P#chr$N.}
    for T in htn dm lip tchl; do echo "$T $N $P" >> ${PROJECT_ROOT}/work/run_wg4/worklist.txt; done
  done
done
echo "worklist: $(wc -l < ${PROJECT_ROOT}/work/run_wg4/worklist.txt) runs"
xargs -P "$POOL" -L1 bash -c 'run_one "$@"' _ < ${PROJECT_ROOT}/work/run_wg4/worklist.txt

n=$(ls $OUT/*.done 2>/dev/null | wc -l)
echo "completed: $n/$(( NCH * 4 ))"
[ "$n" -eq "$(( NCH * 4 ))" ] || { echo "GATE FAIL: $n/$(( NCH * 4 ))"; exit 6; }
echo "WG_V4_COMPLETE"
