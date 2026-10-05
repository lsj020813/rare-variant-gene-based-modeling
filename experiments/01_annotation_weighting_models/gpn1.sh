#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
N="${1:?}"; R=${PROJECT_ROOT}/work/ref15
M=$R/annot/map38/chr$N.map.tsv
S=$R/gpnmsa/scores.tsv.bgz
OUT=$R/annot/gpn; mkdir -p $OUT; O=$OUT/chr$N.gpn.tsv
TB=$(ls tabix 2>/dev/null | head -1)
[ -n "$TB" ] || { echo GATE FAIL: no tabix; exit 3; }
MROWS=$(wc -l < $M)
awk -F'\t' '{split($1,a,":"); print a[1]"\t"a[2]-1"\t"a[2]}' $M | sort -k1,1 -k2,2n -u > $O.bed
$TB -R $O.bed $S | LC_ALL=C awk -F'\t' -v M=$M -v OUT=$O '
BEGIN{ while((getline l < M)>0){ split(l,a,"\t"); K[a[1]]=a[2] }
       print "key37\tgpn_msa" > OUT }
{ k=$1":"$2":"$3":"$4; if(k in K){ hit++; print K[k]"\t"$5 >> OUT } }
END{ printf "GPNSTATS\t%d\n", hit > "/dev/stderr" }'
L=$(( $(wc -l < $O) - 1 ))
[ $L -gt 0 ] || { echo GATE FAIL: 0 rows; exit 5; }
[ $L -le $MROWS ] || { echo GATE FAIL: over map; exit 6; }
awk -v a=$L -v b=$MROWS 'BEGIN{printf "[chr%s] gpn %d / map %d = %.1f%%\n", "'$N'", a, b, a/b*100}'
rm -f $O.bed; echo "ok $L" > $O.done; echo GPN_DONE
