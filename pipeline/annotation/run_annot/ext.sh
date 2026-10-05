#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
N="${1:?chr number required}"
R=${PROJECT_ROOT}/work/ref
F=$R/favor/chr$N.tar.gz
M=$R/annot/map38/chr$N.map.tsv
OUT=$R/annot/extract; mkdir -p "$OUT"
O=$OUT/chr$N.annot.tsv
[ -s "$F" ] || { echo "[chr$N] GATE FAIL: no FAVOR $F"; exit 3; }
[ -s "$M" ] || { echo "[chr$N] GATE FAIL: no map $M"; exit 3; }
MROWS=$(wc -l < "$M")
echo "[chr$N] map $MROWS rows"
S=$(date +%s)
tar xzOf "$F" 2>/dev/null | LC_ALL=C gawk -F, -v M="$M" -v CH="$N" -v OUT="$O" '
BEGIN{ while((getline l < M) > 0){ split(l, a, "\t"); K[a[1]] = a[2] }
       print "key37\tcons\tepi_active\tepi_repr\tepi_trans\ttf\tcage_prom\tgenehancer\tlinsight" > OUT }
NR == 1 { next }
{ n++; off = NF - 35
  if (off < 0) { short++; next }
  if (off > maxoff) maxoff = off
  k = $1; gsub(/-/, ":", k)
  if (!(k in K)) next
  hit++; gh = $(32 + off); lin = $(33 + off)
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", K[k], $6, $9, $10, $11, $22, $23, gh, lin >> OUT
  if ($6 != "") c_cons++;  if ($9 != "") c_act++;   if ($22 != "") c_tf++
  if ($23 != "") c_cage++; if (gh != "") c_gh++;    if (lin != "") c_lin++ }
END{ printf "[chr%s] scanned %d short %d maxoff %d hit %d\n", CH, n, short, maxoff, hit > "/dev/stderr"
     printf "STATS\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", CH,n,short,hit,c_cons,c_act,c_tf,c_cage,c_gh,c_lin > "/dev/stderr" }'
RC=$?
[ "$RC" -eq 0 ] || { echo "[chr$N] GATE FAIL: pipe rc=$RC"; exit 4; }
LINES=$(( $(wc -l < "$O") - 1 ))
echo "[chr$N] out $LINES rows (map $MROWS) elapsed $(( $(date +%s) - S ))s"
[ "$LINES" -gt 0 ]        || { echo "[chr$N] GATE FAIL: 0 rows"; exit 5; }
[ "$LINES" -le "$MROWS" ] || { echo "[chr$N] GATE FAIL: $LINES > $MROWS"; exit 6; }
awk -v a="$LINES" -v b="$MROWS" 'BEGIN{printf "[chr'"$N"'] join %.1f%%\n", a/b*100}'
echo "ok $LINES" > "$O.done"
echo "[chr$N] EXTRACT_DONE"
