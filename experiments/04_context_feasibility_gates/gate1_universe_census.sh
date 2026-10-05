#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
B=bcftools
W=${PROJECT_ROOT}/work/gate1; mkdir -p $W/out $W/logs
echo "START $(date -Is) pid=$$ host=$(hostname)"
for N in 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1; do
  V=${PROJECT_ROOT}/work/ref/orig_index/chr$N.vcf.gz
  [ -s "$V" ] || { echo "SKIP chr$N"; continue; }
  O=$W/out/universe_chr$N.tsv
  [ -s "$O" ] && { echo "HAVE chr$N"; continue; }
  nice -n 19 ionice -c3 $B query -f '%INFO/MAF\t%INFO/R2\t%INFO/TYPED\t%INFO/IMPUTED\t%REF\t%ALT\n' $V \
  | awk -F'\t' '
    { m=$1+0; r=$2+0;
      b=(m<0.001)?"B1_lt0.1":(m<0.01)?"B2_0.1to1":(m<0.05)?"B3_1to5":"B4_ge5";
      rb=(r<0.3)?"R2_lt0.3":(r<0.6)?"R2_0.3to0.6":(r<0.8)?"R2_0.6to0.8":(r<0.9)?"R2_0.8to0.9":"R2_ge0.9";
      t=($3=="1"||$3==".1")?"TYPED":"IMPUTED";
      bi=(length($5)==1 && length($6)==1)?"SNV":"INDEL";
      multi=(index($6,",")>0)?"MULTI":"BI";
      key=b"\t"rb"\t"t"\t"bi"\t"multi; n[key]++; s[key]+=r; N0++ }
    END{ for(k in n) printf "%s\t%d\t%.4f\n", k, n[k], s[k]/n[k]; printf "TOTAL\t.\t.\t.\t.\t%d\t.\n", N0 }' > $O.tmp && mv $O.tmp $O
  echo "DONE chr$N $(date -Is) $(wc -l < $O) rows"
done
echo "END $(date -Is)"
echo "gate1 universe census complete $(date -Is)" > $W/out/universe.done
