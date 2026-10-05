#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
N="${1:?chr number required}"
R=${PROJECT_ROOT}/work/ref
F=$R/favor/chr$N.tar.gz
M=$R/annot/map38_trA/chr$N.map.tsv
OUT=$R/annot/extract_trA; mkdir -p "$OUT"
O=$OUT/chr$N.annot.tsv; OC=$OUT/chr$N.cadd.tsv
[ -s "$F" ] || { echo "[chr$N] GATE FAIL: no FAVOR $F"; exit 3; }
[ -s "$M.done" ] || { echo "[chr$N] GATE FAIL: no map $M"; exit 3; }
[ -s "$O.done" ] && [ -s "$OC.done" ] && { echo "[chr$N] 이미 완료"; exit 0; }
MROWS=$(wc -l < "$M")
echo "[chr$N] map $MROWS rows"
S=$(date +%s)
tar xzOf "$F" 2>/dev/null | LC_ALL=C gawk -F, -v M="$M" -v CH="$N" -v OUT="$O.tmp" -v OUTC="$OC.tmp" '
BEGIN{ while((getline l < M) > 0){ split(l, a, "\t"); K[a[1]] = a[2] }
       print "key37\tcons\tepi_active\tepi_repr\tepi_trans\ttf\tcage_prom\tgenehancer\tlinsight" > OUT }
NR == 1 { printf "HEADER\tNF=%d\t1=%s\t6=%s\t9=%s\t10=%s\t11=%s\t22=%s\t23=%s\t29=%s\t32=%s\t33=%s\t34=%s\t35=%s\n", NF,$1,$6,$9,$10,$11,$22,$23,$29,$32,$33,$34,$35 > "/dev/stderr"; next }
{ n++
  if (NF < 35) { short++; next }
  rd = 0
  if ($NF ~ /"$/) { j = NF; while (j > 1 && $j !~ /^"/) j--; rd = NF - j; if (rd > 0) rdq++ }
  off = NF - 35 - rd
  if (off < 0) { short++; next }
  if (off > maxoff) maxoff = off
  k = $1; gsub(/-/, ":", k)
  if (!(k in K)) next
  hit++; gh = $(32 + off); lin = $(33 + off)
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", K[k], $6, $9, $10, $11, $22, $23, gh, lin >> OUT
  cd = $(NF - 1 - rd); printf "%s\t%s\n", K[k], cd >> OUTC; if (cd != "") c_cadd++
  if ($6 != "") c_cons++;  if ($9 != "") c_act++;   if ($22 != "") c_tf++
  if ($23 != "") c_cage++; if (gh != "") c_gh++;    if (lin != "") c_lin++ }
END{ printf "[chr%s] scanned %d short %d maxoff %d rdhs_quoted_comma_rows %d hit %d\n", CH, n, short, maxoff, rdq, hit > "/dev/stderr"
     printf "STATS\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", CH,n,short,hit,c_cons,c_act,c_tf,c_cage,c_gh,c_lin,c_cadd > "/dev/stderr" }'
RC=$?
[ "$RC" -eq 0 ] || { echo "[chr$N] GATE FAIL: pipe rc=$RC"; exit 4; }
LINES=$(( $(wc -l < "$O.tmp") - 1 ))
LC=$(wc -l < "$OC.tmp")
echo "[chr$N] out $LINES rows cadd $LC rows (map $MROWS) elapsed $(( $(date +%s) - S ))s"
[ "$LINES" -gt 0 ]        || { echo "[chr$N] GATE FAIL: 0 rows"; exit 5; }
[ "$LINES" -le "$MROWS" ] || { echo "[chr$N] GATE FAIL: $LINES > $MROWS"; exit 6; }
[ "$LC" -eq "$LINES" ]    || { echo "[chr$N] GATE FAIL: cadd rows $LC != annot rows $LINES"; exit 6; }
NUM=$(awk -F'\t' '$2!="" && $2 !~ /^[0-9.]+$/' "$OC.tmp" | wc -l)
[ "$NUM" -eq 0 ] || { echo "[chr$N] GATE FAIL: cadd 비수치 값 $NUM 행"; exit 7; }
awk -v a="$LINES" -v b="$MROWS" 'BEGIN{printf "[chr'"$N"'] join %.1f%%\n", a/b*100}'
mv "$O.tmp" "$O"; mv "$OC.tmp" "$OC"
echo "ok $LINES" > "$O.done"; echo "ok $LC" > "$OC.done"
echo "[chr$N] EXTRACT_DONE"
