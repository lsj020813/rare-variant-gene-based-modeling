#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -u
BCF=bcftools
R=${PROJECT_ROOT}/work/ref; OUT=$R/band15_vcf; mkdir -p $OUT
one() {
  N=$1; S=$R/orig_index/chr$N.vcf.gz; O=$OUT/chr$N.band15.vcf.gz
  [ -s $O.done ] && { echo "[chr$N] skip (done)"; return; }
  [ -s $S ] || { echo "[chr$N] GATE FAIL: source missing"; return; }
  echo "[chr$N] start $(date '+%H:%M')"
  $BCF view -i 'INFO/MAF>0.01 && INFO/MAF<=0.05' -Oz -o $O.tmp $S && mv $O.tmp $O && $BCF index -c $O || { echo "[chr$N] GATE FAIL: view/index"; rm -f $O.tmp; return; }
  n=$($BCF index -n $O)
  [ "$n" -gt 0 ] || { echo "[chr$N] GATE FAIL: 0 records"; return; }
  echo "ok $n" > $O.done; echo "[chr$N] done $n records $(date '+%H:%M')"
}
export -f one; export BCF R OUT
printf '%s\n' 15 16 17 18 19 20 21 22 | xargs -P 8 -I{} bash -c 'one {}'
echo "completed: $(ls $OUT/*.done 2>/dev/null | wc -l)/22"
echo BAND15_EXTRACT_DONE
