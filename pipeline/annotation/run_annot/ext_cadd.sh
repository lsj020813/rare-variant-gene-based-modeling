#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
N=$1; F=${PROJECT_ROOT}/work/ref/favor/chr$N.tar.gz; M=${PROJECT_ROOT}/work/ref/annot/map38/chr$N.map.tsv
O=${PROJECT_ROOT}/work/ref/annot/extract/chr$N.cadd.tsv
[ -s "$F" ] && [ -s "$M" ] || { echo "[chr$N] GATE FAIL: 입력 없음"; exit 2; }
[ -s "$O.done" ] && { echo "[chr$N] 이미 완료"; exit 0; }
MROWS=$(wc -l < "$M"); S=$(date +%s)
tar xzOf "$F" 2>/dev/null | LC_ALL=C awk -F, -v M="$M" -v OUT="$O.tmp" '
  BEGIN{while((getline l < M)>0){split(l,a,"\t"); K[a[1]]=a[2]}}
  NR>1{ if(NF<35){bad++; next} k=$1; gsub(/-/,":",k);
        if(k in K){ hit++; if($(NF-1)!="") nn++; printf "%s\t%s\n", K[k], $(NF-1) > OUT } }
  END{printf "[chr'$N'] 행 %d | NF<35 %d | 히트 %d | cadd 비결측 %d\n", NR-1, bad, hit, nn}'
LINES=$(wc -l < "$O.tmp" 2>/dev/null || echo 0)
[ "$LINES" -gt 0 ] || { echo "[chr$N] GATE FAIL: 출력 0행"; exit 5; }
[ "$LINES" -le "$MROWS" ] || { echo "[chr$N] GATE FAIL: $LINES > 대응표 $MROWS"; exit 6; }
NUM=$(awk -F'\t' '$2!="" && $2 !~ /^[0-9.]+$/' "$O.tmp" | head -3 | wc -l)
[ "$NUM" -eq 0 ] || { echo "[chr$N] GATE FAIL: cadd 비수치 값 존재"; awk -F'\t' '$2!="" && $2 !~ /^[0-9.]+$/' "$O.tmp" | head -3; exit 7; }
mv "$O.tmp" "$O"; echo "ok $LINES" > "$O.done"
echo "[chr$N] CADD_DONE $LINES 행 $(( $(date +%s) - S ))초"
