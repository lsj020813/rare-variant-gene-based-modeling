#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
exec 9>${PROJECT_ROOT}/work/run_bwg2/.lock
flock -n 9 || { echo "already running"; exit 0; }

R=${PROJECT_ROOT}/work/ref
S1=$R/saige_step1_v4
GP=$R/groupfiles_bwg
BAND=$R/band_vcf
OUT=$R/saige_step2_bwg
mkdir -p "$OUT"
export UDOCKER_DIR=${CONTAINER_STATE_DIR:?Set CONTAINER_STATE_DIR}
UD="PROOT_NO_SECCOMP=1 ${UDOCKER_BIN:-udocker} run --volume=/data:/data ${SAIGE_IMAGE:?Set SAIGE_IMAGE}"
POOL=${POOL:-16}

run_one() {
  local T=$1 N=$2
  local G=$GP/chr$N.B_3kb_re2g.txt
  local O=$OUT/$T.chr$N
  [ -s "$O.done" ] && { echo "  skip $T chr$N"; return 0; }
  [ -s "$G" ] || { echo "  MISS $T chr$N (no groupfile)"; return 1; }
  eval $UD step2_SPAtests.R \
    --vcfFile="$BAND/chr$N.band.vcf.gz" --vcfFileIndex="$BAND/chr$N.band.vcf.gz.csi" \
    --vcfField="DS" --chrom="$N" --AlleleOrder=ref-first \
    --minMAF=0 --minMAC=0.5 --LOCO=FALSE --is_fastTest=TRUE \
    --GMMATmodelFile="$S1/${T}_v4.rda" --varianceRatioFile="$S1/${T}_v4.varianceRatio.txt" \
    --groupFile="$G" --annotation_in_groupTest="all" --maxMAF_in_groupTest=0.01 \
    --SAIGEOutputFile="$O" > "$O.log" 2>&1
  local rc=$?
  local exp=$(awk '$2=="var"{print $1}' "$G" | sort -u | wc -l)
  local nrow=0; [ -s "$O" ] && nrow=$(( $(wc -l < "$O") - 1 ))
  if [ "$rc" -ne 0 ] || [ "$nrow" -le 0 ]; then
    echo "  FAIL $T chr$N rc=$rc rows=$nrow/$exp"; tail -2 "$O.log"; return 1
  fi
  echo ok > "$O.done"
  [ "$nrow" -eq "$exp" ] && echo "  OK $T chr$N rows=$nrow" || echo "  PARTIAL $T chr$N $nrow/$exp"
}
export -f run_one; export GP BAND OUT S1 UD R

: > $OUT/joblist.txt
for N in $(seq 22 -1 1); do for T in tchl htn dm lip; do echo "$T $N" >> $OUT/joblist.txt; done; done
cat $OUT/joblist.txt | xargs -P "$POOL" -n2 bash -c 'run_one "$@"' _

nd=$(ls $OUT/*.done 2>/dev/null | wc -l)
echo "completed: $nd / 88"
echo "BWG_ASSOC_COMPLETE"
