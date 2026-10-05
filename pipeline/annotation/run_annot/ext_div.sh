#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
N=$1
F=${PROJECT_ROOT}/work/ref/favor/chr$N.tar.gz
M=${PROJECT_ROOT}/work/ref/annot/map38/chr$N.map.tsv
O=${PROJECT_ROOT}/work/ref/annot/extract/chr$N.div.tsv
[ -s "$F" ] || { echo "[chr$N] GATE FAIL: no FAVOR $F"; exit 3; }
[ -s "$M" ] || { echo "[chr$N] GATE FAIL: no map $M"; exit 3; }
rm -f "$O" "$O.done"
S=$(date +%s)
tar xzOf "$F" 2>/dev/null | LC_ALL=C gawk -F, -v M="$M" -v OUT="$O" '
BEGIN{ while((getline l < M)>0){ split(l,a,"\t"); K[a[1]]=a[2] }
       print "key\tdiv1\tdiv2\tdiv3\tmutdens\tcons_v2" > OUT }
NR>1{ k=$1; gsub(/-/,":",k); if(!(k in K)) next
      if($1 !~ /^[0-9XY]+-[0-9]+-[ACGT]+-[ACGT]+$/) { bad++; next }
      printf "%s\t%s\t%s\t%s\t%s\t%s\n", K[k], $12, $13, $14, $17, $7 >> OUT
      n++; if($12!="") c1++; if($14!="") c3++; if($17!="") cm++ }
END{ printf "[extract] rows %d | div1 nonnull %d | div3 nonnull %d | mutdens nonnull %d | bad_key %d\n", n, c1, c3, cm, bad > "/dev/stderr" }'
RC=${PIPESTATUS[1]}
[ "$RC" -eq 0 ] || { echo "[chr$N] GATE FAIL: awk rc=$RC"; exit 4; }
LINES=$(( $(wc -l < "$O") - 1 )); MROWS=$(wc -l < "$M")
[ "$LINES" -gt 0 ]        || { echo "[chr$N] GATE FAIL: 0 rows"; exit 5; }
[ "$LINES" -le "$MROWS" ] || { echo "[chr$N] GATE FAIL: $LINES > $MROWS"; exit 6; }
awk -v a="$LINES" -v b="$MROWS" 'BEGIN{printf "[chr'"$N"'] join %.1f%%\n", a/b*100}'
echo "ok $LINES" > "$O.done"
echo "[chr$N] DIV_EXTRACT_DONE $(( $(date +%s) - S ))s"
