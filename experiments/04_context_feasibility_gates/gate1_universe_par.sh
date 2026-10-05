#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
N=$1; SRC=$2; TAG=$3
B=bcftools
W=${PROJECT_ROOT}/work/gate1
O=$W/out/universe_${TAG}_chr$N.tsv
[ -s "$O" ] && { echo "HAVE $TAG chr$N"; exit 0; }
V=$SRC
[ -s "$V" ] || { echo "MISSING $V"; exit 1; }
nice -n 19 ionice -c3 $B query -f '%INFO/MAF\t%INFO/R2\t%INFO/TYPED\t%REF\t%ALT\n' "$V" \
| awk -F'\t' -v chr=$N -v tag=$TAG '
  { m=$1+0; r=$2+0;
    b=(m<0.001)?"B1_lt0.1pct":(m<0.01)?"B2_0.1to1pct":(m<0.05)?"B3_1to5pct":"B4_ge5pct";
    rb=(r<0.3)?"R2lt0.3":(r<0.6)?"R2_0.3to0.6":(r<0.8)?"R2_0.6to0.8":(r<0.9)?"R2_0.8to0.9":"R2ge0.9";
    t=($3=="1")?"TYPED":"IMPUTED";
    vt=(length($4)==1 && length($5)==1 && $5!~/,/)?"SNV":"INDELorMNP";
    al=(index($5,",")>0)?"MULTI":"BI";
    n[b"\t"rb"\t"t"\t"vt"\t"al]++; tot++ }
  END{ for(k in n) printf "%s\t%s\t%s\t%d\n", tag, chr, k, n[k]; printf "%s\t%s\tTOTAL\t.\t.\t.\t.\t%d\n", tag, chr, tot }' > $O.tmp && mv $O.tmp $O
echo "DONE $TAG chr$N $(date -Is)"
