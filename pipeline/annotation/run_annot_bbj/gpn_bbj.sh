#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
N="${1:?}"; R=${PROJECT_ROOT}/work/ref
M=$R/annot/map38_bbj/chr$N.map.tsv
S=$R/gpnmsa/scores.tsv.bgz
OUT=$R/annot/gpn_bbj; mkdir -p $OUT; O=$OUT/chr$N.gpn.tsv
[ -s "$O.done" ] && { echo "[chr$N] 이미 완료"; exit 0; }
[ -s "$M.done" ] || { echo GATE FAIL: no map; exit 3; }
TB=$(ls tabix 2>/dev/null | head -1)
[ -n "$TB" ] || { echo GATE FAIL: no tabix; exit 3; }
MROWS=$(wc -l < $M)
awk -F'\t' '{split($1,a,":"); print a[1]"\t"a[2]-1"\t"a[2]}' $M | sort -k1,1 -k2,2n -u > $O.bed
$TB -R $O.bed $S | LC_ALL=C awk -F'\t' -v M=$M -v OUT=$O.tmp '
BEGIN{ while((getline l < M)>0){ split(l,a,"\t"); K[a[1]]=a[2] }
       print "key37\tgpn_msa" > OUT }
{ k=$1":"$2":"$3":"$4; if(k in K){ hit++; print K[k]"\t"$5 >> OUT } }
END{ printf "GPNSTATS\t%d\n", hit > "/dev/stderr" }'
RC=$?; [ $RC -eq 0 ] || { echo GATE FAIL: pipe rc=$RC; exit 4; }
L=$(( $(wc -l < $O.tmp) - 1 ))
[ $L -gt 0 ] || { echo GATE FAIL: 0 rows; exit 5; }
[ $L -le $MROWS ] || { echo GATE FAIL: over map; exit 6; }
awk -v a=$L -v b=$MROWS 'BEGIN{printf "[chr%s] gpn %d / map %d = %.1f%%\n", "'$N'", a, b, a/b*100}'
mv $O.tmp $O; rm -f $O.bed; echo "ok $L" > $O.done; echo GPN_DONE
