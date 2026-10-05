#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
: "${SAIGE_IMAGE:?Set SAIGE_IMAGE}"
: "${UDOCKER_BIN:?Set UDOCKER_BIN}"
set -uo pipefail
N=$1; R=${PROJECT_ROOT}/work/ref; V=$R/common05/chr$N.maf05.vcf.gz; O=$R/gwas05; mkdir -p $O
UD=${UDOCKER_BIN}
for i in $(seq 1 1440); do [ -s "$V.done" ] && break; sleep 30; done
[ -s "$V.done" ] || { echo "[chr$N] GATE FAIL: 추출 마커 대기 초과"; exit 2; }
[ -s "$V.csi" ] || { echo "[chr$N] GATE FAIL: csi 없음"; exit 3; }
for T in tchl htn dm lip; do
  OUT=$O/$T.chr$N.txt
  [ -s "$OUT.done" ] && { echo "[chr$N $T] 이미 완료"; continue; }
  S=$(date +%s)
  $UD run --volume=/data:/data ${SAIGE_IMAGE} \
    step2_SPAtests.R \
      --vcfFile="$V" --vcfFileIndex="$V.csi" --vcfField=DS \
      --chrom=$N --AlleleOrder=ref-first --minMAF=0.05 --minMAC=20 \
      --GMMATmodelFile="$R/saige_step1_v4/${T}_v4.rda" \
      --varianceRatioFile="$R/saige_step1_v4/${T}_v4.varianceRatio.txt" \
      --LOCO=FALSE --SAIGEOutputFile="$OUT" > "$OUT.log" 2>&1
  RC=$?
  [ -s "$OUT" ] || { echo "[chr$N $T] GATE FAIL: 출력 없음 rc=$RC"; tail -3 "$OUT.log"; exit 4; }
  head -1 "$OUT" | grep -q $'^CHR\tPOS\tMarkerID' || { echo "[chr$N $T] GATE FAIL: 헤더 불일치"; exit 5; }
  NR_=$(( $(wc -l < "$OUT") - 1 )); [ "$NR_" -gt 1000 ] || { echo "[chr$N $T] GATE FAIL: 행 $NR_"; exit 6; }
  NREC=$(cut -d' ' -f2 "$V.done"); PCT=$(awk -v a=$NR_ -v b=$NREC 'BEGIN{printf "%.1f", a/b*100}')
  echo "ok $NR_" > "$OUT.done"
  echo "[chr$N $T] GWAS_DONE 행 $NR_ (변이 $NREC 의 ${PCT}%) $(( $(date +%s) - S ))s"
done
echo "[chr$N] CHR_DONE"
