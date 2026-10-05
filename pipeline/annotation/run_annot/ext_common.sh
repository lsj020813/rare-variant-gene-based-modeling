#!/usr/bin/env bash
: "${GENOTYPE_DIR:?Set GENOTYPE_DIR}"
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
N=$1; RAW=${GENOTYPE_DIR}/chr$N.vcf.gz
D=${PROJECT_ROOT}/work/ref/common05; mkdir -p $D; O=$D/chr$N.maf05.vcf.gz
BCF=bcftools
[ -s "$RAW" ] || { echo "[chr$N] GATE FAIL: 원본 없음"; exit 2; }
[ -s "$O.done" ] && { echo "[chr$N] 이미 완료"; exit 0; }
S=$(date +%s)
$BCF view -i 'INFO/MAF>=0.05' --threads 6 -Ou "$RAW" | $BCF annotate -x FORMAT/GT,FORMAT/GP --threads 6 -Oz -o "$O.tmp" || { echo "[chr$N] GATE FAIL: bcftools 파이프 종료코드 $?"; rm -f "$O.tmp"; exit 3; }
FMT=$($BCF view -H "$O.tmp" 2>/dev/null | head -1 | cut -f9)
[ "$FMT" = "DS" ] || { echo "[chr$N] GATE FAIL: FORMAT='$FMT' != DS"; exit 8; }
mv "$O.tmp" "$O"; $BCF index -f -c --threads 2 "$O" || { echo "[chr$N] GATE FAIL: index"; exit 4; }
NREC=$($BCF index -n "$O"); [ "$NREC" -gt 1000 ] || { echo "[chr$N] GATE FAIL: 레코드 $NREC"; exit 5; }
MINMAF=$($BCF query -f '%INFO/MAF\n' "$O" | head -20000 | sort -g | head -1)
awk -v m="$MINMAF" 'BEGIN{exit !(m>=0.05)}' || { echo "[chr$N] GATE FAIL: MAF 최소 $MINMAF < 0.05"; exit 6; }
echo "ok $NREC" > "$O.done"
echo "[chr$N] COMMON_DONE $NREC 변이 $(( $(date +%s) - S ))초 $(du -h $O | cut -f1)"
