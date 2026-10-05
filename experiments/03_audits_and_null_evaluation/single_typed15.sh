#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
exec 9>${PROJECT_ROOT}/work/run_meth/.lock_single_typed
flock -n 9 || { echo "already running"; exit 0; }
R=${PROJECT_ROOT}/work/ref; R15=${PROJECT_ROOT}/work/ref15
S1=$R/saige_step1_v4; BAND=$R15/band_vcf
OUT=$R15/saige_single_typed; RNG=$OUT/ranges; mkdir -p "$OUT" "$RNG"
BCF=${BCFTOOLS:-bcftools}
export UDOCKER_DIR=${CONTAINER_STATE_DIR:?Set CONTAINER_STATE_DIR}
UD="PROOT_NO_SECCOMP=1 ${UDOCKER_BIN:-udocker} run --volume=/data:/data ${SAIGE_IMAGE:?Set SAIGE_IMAGE}"
POOL=${POOL:-8}
echo "START $(date -Is) pid=$$ host=$(hostname)"

build_rng() { local N=$1; local F=$RNG/chr$N.typed.ranges
  [ -s "$F.done" ] && return 0
  $BCF query -i 'INFO/TYPED=1' -f '%CHROM\t%POS\t%POS\n' "$BAND/chr$N.band.vcf.gz" > "$F" && [ -s "$F" ] && echo "$(wc -l < $F)" > "$F.done" && echo "  ranges chr$N n=$(cat $F.done)"; }
export -f build_rng; export BCF BAND RNG
seq 22 -1 1 | xargs -P 4 -n1 bash -c 'build_rng "$@"' _
echo "RANGES_DONE $(date -Is) total=$(cat $RNG/*.done | paste -sd+ | bc)"

run_one() { local T=$1 N=$2; local O=$OUT/$T.chr$N; local F=$RNG/chr$N.typed.ranges
  [ -s "$O.done" ] && { echo "  skip $T chr$N"; return 0; }
  [ -s "$F" ] || { echo "  MISS $T chr$N (no ranges)"; return 1; }
  eval $UD step2_SPAtests.R \
    --vcfFile="$BAND/chr$N.band.vcf.gz" --vcfFileIndex="$BAND/chr$N.band.vcf.gz.csi" \
    --vcfField="DS" --chrom="$N" --AlleleOrder=ref-first \
    --rangestoIncludeFile="$F" \
    --minMAF=0 --minMAC=20 --LOCO=FALSE --is_fastTest=TRUE \
    --GMMATmodelFile="$S1/${T}_v4.rda" --varianceRatioFile="$S1/${T}_v4.varianceRatio.txt" \
    --SAIGEOutputFile="$O" > "$O.log" 2>&1
  local rc=$?; local nrow=0; [ -s "$O" ] && nrow=$(( $(wc -l < "$O") - 1 ))
  if [ "$rc" -ne 0 ] || [ "$nrow" -le 0 ]; then echo "  FAIL $T chr$N rc=$rc rows=$nrow"; tail -2 "$O.log"; return 1; fi
  echo "rows=$nrow $(date -Is)" > "$O.done"; echo "  OK $T chr$N rows=$nrow exp=$(cat $F.done)"; }
export -f run_one; export OUT RNG S1 UD BAND
: > $OUT/joblist.txt
for N in $(seq 22 -1 1); do for T in tchl htn dm lip; do echo "$T $N" >> $OUT/joblist.txt; done; done
cat $OUT/joblist.txt | xargs -P "$POOL" -n2 bash -c 'run_one "$@"' _
nd=$(ls $OUT/*.done 2>/dev/null | grep -v ranges | wc -l)
echo "completed: $nd / 88  $(date -Is)"
[ "$nd" -eq 88 ] && echo "single_typed15 88/88 $(date -Is)" > $OUT/ALL.done
echo "END $(date -Is)"
