#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -u
N=$1; BCF=bcftools
O=${PROJECT_ROOT}/work/ref/common05/chr$N.maf05.vcf.gz
[ -s "$O.done" ] && { echo "[chr$N] 이미 완료"; exit 0; }
[ -s "$O.tmp" ] || { echo "[chr$N] tmp 없음"; exit 2; }
pgrep -f "ext_common.sh $N\$" >/dev/null && { echo "[chr$N] 워커 실행 중 — 대기"; exit 3; }
FMT=$($BCF view -H "$O.tmp" 2>/dev/null | head -1 | cut -f9)
[ "$FMT" = "DS" ] || { echo "[chr$N] GATE FAIL: FORMAT='$FMT' != DS"; exit 8; }
mv "$O.tmp" "$O"; $BCF index -f -c --threads 2 "$O" || { echo "[chr$N] GATE FAIL: index"; exit 4; }
NREC=$($BCF index -n "$O"); [ "$NREC" -gt 1000 ] || { echo "[chr$N] GATE FAIL: 레코드 $NREC"; exit 5; }
MINMAF=$($BCF query -f '%INFO/MAF\n' "$O" | sort -g | head -1)
awk -v m="$MINMAF" 'BEGIN{exit !(m>=0.05)}' || { echo "[chr$N] GATE FAIL: MAF 최소 $MINMAF < 0.05"; exit 6; }
echo "ok $NREC" > "$O.done"
echo "[chr$N] FINISH_DONE 레코드 $NREC 최소MAF $MINMAF $(du -h $O | cut -f1)"
