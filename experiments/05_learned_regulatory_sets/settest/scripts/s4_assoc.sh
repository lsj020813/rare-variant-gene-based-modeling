#!/usr/bin/env bash
: "${CONTAINER_STATE_DIR:?Set CONTAINER_STATE_DIR}"
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
: "${SAIGE_IMAGE:?Set SAIGE_IMAGE}"
: "${UDOCKER_BIN:?Set UDOCKER_BIN}"
set -o pipefail
BASE=${PROJECT_ROOT}/work/fset/settest
SUB=$BASE/sub; GP=$BASE/groupfiles; OUT=$BASE/assoc; LOG=$BASE/logs
mkdir -p "$OUT" "$LOG"
export TMPDIR=${PROJECT_ROOT}/work/tmp TMP=${PROJECT_ROOT}/work/tmp TEMP=${PROJECT_ROOT}/work/tmp
export UDOCKER_DIR=${CONTAINER_STATE_DIR}
S1=${PROJECT_ROOT}/work/ref/saige_step1_v4
UD="PROOT_NO_SECCOMP=1 ${UDOCKER_BIN} run --volume=/data:/data ${SAIGE_IMAGE}"

GTAG=${GTAG:-set}
MAXMAF=${MAXMAF:-0.5,0.01}
TRAITS=${TRAITS:-"tchl htn dm lip"}
CHRS=${CHRS:-"$(seq 22 -1 1)"}
POOL=${POOL:-6}

run_one() {
  local T=$1 N=$2
  local G=$GP/chr$N.$GTAG.txt
  local O=$OUT/$GTAG.$T.chr$N
  [ -s "$O.done" ] && { echo "  skip $T chr$N"; return 0; }
  [ -s "$G" ] || { echo "  MISS $T chr$N (no groupfile)"; return 1; }
  local t0=$(date +%s)
  eval nice -n 19 ionice -c3 env $UD step2_SPAtests.R \
    --vcfFile="$SUB/chr$N.sub.vcf.gz" --vcfFileIndex="$SUB/chr$N.sub.vcf.gz.csi" \
    --vcfField="DS" --chrom="$N" --AlleleOrder=ref-first \
    --minMAF=0 --minMAC=0.5 --LOCO=FALSE --is_fastTest=TRUE \
    --GMMATmodelFile="$S1/${T}_v4.rda" --varianceRatioFile="$S1/${T}_v4.varianceRatio.txt" \
    --groupFile="$G" --annotation_in_groupTest="all" --maxMAF_in_groupTest="$MAXMAF" \
    --SAIGEOutputFile="$O" > "$O.log" 2>&1
  local rc=$?
  local exp=$(awk '$2=="var"{print $1}' "$G" | sort -u | wc -l)
  local nrow=0; [ -s "$O" ] && nrow=$(( $(wc -l < "$O") - 1 ))
  local t1=$(date +%s)
  if [ "$rc" -ne 0 ] || [ "$nrow" -le 0 ]; then
    echo "  FAIL $T chr$N rc=$rc rows=$nrow groups=$exp secs=$((t1-t0))"; tail -3 "$O.log"; return 1
  fi
  echo "rows=$nrow groups=$exp secs=$((t1-t0))" > "$O.done"
  echo "  OK $T chr$N rows=$nrow groups=$exp secs=$((t1-t0))"
}
export -f run_one; export GP SUB OUT S1 UD GTAG MAXMAF

JL=$OUT/joblist.$GTAG.txt; : > "$JL"
for N in $CHRS; do for T in $TRAITS; do echo "$T $N" >> "$JL"; done; done
NJ=$(wc -l < "$JL")
cat "$JL" | xargs -P "$POOL" -n2 bash -c 'run_one "$@"' _
nd=$(ls $OUT/$GTAG.*.done 2>/dev/null | wc -l)
echo "completed: $nd / $NJ"
[ "$nd" -eq "$NJ" ] && echo "ASSOC_COMPLETE $GTAG" > $BASE/assoc.$GTAG.done
