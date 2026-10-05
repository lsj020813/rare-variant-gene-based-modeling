#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
B=${BCFTOOLS:-bcftools}
W=${PROJECT_ROOT}/work/run_meth/out; mkdir -p $W
echo "START $(date -Is) pid=$$"
for band in ref ref15; do
 for N in $(seq 1 22); do
  V=${PROJECT_ROOT}/work/$band/band_vcf/chr$N.band.vcf.gz
  [ -s "$V" ] || { echo "$band chr$N MISSING"; continue; }
  nice -n 19 ionice -c3 $B query -f '%INFO/MAF\t%INFO/R2\t%INFO/TYPED\n' $V | awk -F'\t' -v band=$band -v chr=$N '
   { m=$1+0; b=(m<0.001)?"a<0.1%":(m<0.005)?"b0.1-0.5%":(m<0.01)?"c0.5-1%":(m<0.02)?"d1-2%":(m<0.05)?"e2-5%":"f>=5%";
     typed=($3!=".")?"TYPED":"IMPUTED"; k=band"\t"chr"\t"b"\t"typed; n[k]++; r2[k]+=$2 }
   END{ for(k in n) printf "%s\t%d\t%.4f\n", k, n[k], r2[k]/n[k] }' | sort
  echo "DONE $band chr$N $(date -Is)"
 done
done
echo "END $(date -Is)"
