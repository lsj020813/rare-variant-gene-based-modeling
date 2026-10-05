#!/usr/bin/env bash
: "${CONTAINER_STATE_DIR:?Set CONTAINER_STATE_DIR}"
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
: "${SAIGE_IMAGE:?Set SAIGE_IMAGE}"
: "${UDOCKER_BIN:?Set UDOCKER_BIN}"
set -o pipefail
TAG=$1; PAR=$2; MAXSTART=$3; MAXRUN=$4
W=${PROJECT_ROOT}/work; B=$W/run_l3b; T=tchl; N=19
R=$W/ref; S1=$R/saige_step1_v4; G=$B/out/chunks; O=$B/out/arm_$TAG
mkdir -p "$O"; rm -f "$B/out/STOP"
export UDOCKER_DIR=${CONTAINER_STATE_DIR}
export TMPDIR=$W/tmp TMP=$W/tmp TEMP=$W/tmp

L1=$(cut -d' ' -f1 /proc/loadavg)
WAITED=0
while [ "$(echo "$L1 > $MAXSTART" | bc -l)" = "1" ]; do
  echo "[$TAG] WAIT $(date '+%H:%M:%S') 1-min load $L1 above start ceiling $MAXSTART (waited ${WAITED}s)"
  sleep 60; WAITED=$((WAITED+60)); L1=$(cut -d' ' -f1 /proc/loadavg)
  if [ "$WAITED" -ge 7200 ]; then echo "[$TAG] HOLD: load stayed above $MAXSTART for 2h"; exit 3; fi
done
[ "$WAITED" -gt 0 ] && echo "[$TAG] WAIT_END $(date '+%H:%M:%S') load $L1 after ${WAITED}s"

( while :; do
    l=$(cut -d' ' -f1 /proc/loadavg)
    if [ "$(echo "$l > $MAXRUN" | bc -l)" = "1" ]; then
      echo "$(date '+%H:%M:%S') load $l exceeded $MAXRUN -- draining" > "$B/out/STOP"; break
    fi
    [ -s "$B/out/$TAG.merged.done" ] && break
    sleep 30
  done ) &
MON=$!

run_one(){
  local P=$1
  [ -s "$B/out/STOP" ] && return 0
  local gf=$G/G_$TAG.$P.txt
  local out=$O/$P
  [ -s $out.done ] && return 0
  PROOT_NO_SECCOMP=1 nice -n 19 ionice -c3 \
    ${UDOCKER_BIN} run --volume=/data:/data ${SAIGE_IMAGE} \
    step2_SPAtests.R \
    --vcfFile=$R/band_vcf/chr$N.band.vcf.gz \
    --vcfFileIndex=$R/band_vcf/chr$N.band.vcf.gz.csi \
    --vcfField=DS --chrom=$N --AlleleOrder=ref-first \
    --minMAF=0 --minMAC=0.5 --LOCO=FALSE --is_fastTest=TRUE --nThreads=1 \
    --GMMATmodelFile=$S1/${T}_v4.rda \
    --varianceRatioFile=$S1/${T}_v4.varianceRatio.txt \
    --groupFile=$gf \
    --annotation_in_groupTest=all --maxMAF_in_groupTest=0.01 \
    --SAIGEOutputFile=$out > $out.log 2>&1
  local rc=$?
  local exp=$(awk '$2=="var"{print $1}' $gf | sort -u | wc -l)
  local nrow=0; [ -s $out ] && nrow=$(( $(wc -l < $out) - 1 ))
  if [ $rc -eq 0 ] && [ "$nrow" -eq "$exp" ]; then echo ok > $out.done
  else echo "FAIL $TAG $P rc=$rc rows=$nrow/$exp"; fi
}
export -f run_one; export G O R S1 T N TAG B UDOCKER_DIR TMPDIR TMP TEMP

NCH=$(ls $G/G_$TAG.part*.txt 2>/dev/null | wc -l)
[ "$NCH" -gt 0 ] || { echo "GATE FAIL: no chunks for $TAG"; kill $MON 2>/dev/null; exit 2; }
echo "[$TAG] start $(date '+%H:%M:%S') chunks=$NCH par=$PAR nThreads=1 load=$L1"
ls $G/G_$TAG.part*.txt | sed 's#.*/G_'"$TAG"'\.##; s#\.txt$##' | \
  xargs -P "$PAR" -n1 bash -c 'run_one "$@"' _
kill $MON 2>/dev/null

ND=$(ls $O/*.done 2>/dev/null | wc -l)
echo "[$TAG] end $(date '+%H:%M:%S') done=$ND/$NCH load=$(cut -d' ' -f1 /proc/loadavg)"
if [ -s "$B/out/STOP" ]; then echo "[$TAG] STOPPED BY LOAD GUARD: $(cat $B/out/STOP)"; exit 4; fi
[ "$ND" -eq "$NCH" ] || { echo "[$TAG] INCOMPLETE"; exit 1; }
head -1 $O/part001 > $B/out/$TAG.merged
for f in $(ls $O/part*[0-9] | sort); do tail -n +2 $f; done >> $B/out/$TAG.merged
NROW=$(( $(wc -l < $B/out/$TAG.merged) - 1 ))
printf '{"status":"complete","tag":"%s","chunks":%d,"rows":%d,"par":%d,"nThreads":1}\n' \
  "$TAG" "$NCH" "$NROW" "$PAR" > $B/out/$TAG.merged.done
echo "[$TAG] merged rows=$NROW"
