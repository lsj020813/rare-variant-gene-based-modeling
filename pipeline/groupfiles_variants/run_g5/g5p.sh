#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
OUT=${PROJECT_ROOT}/work/ref/deductive
mkdir -p "$OUT"
exec 9>"$OUT/.g5p.lock"; flock -n 9 || { echo "another run"; exit 0; }
BCF=${BCFTOOLS:-bcftools}
BED=${BEDTOOLS:-bedtools}
KEYED=${PROJECT_ROOT}/work/ref/lift38_keyed
BAND=${PROJECT_ROOT}/work/ref/band_vcf

one(){
  local N=$1
  [ -e "$OUT/chr${N}.a.done" ] && return 0
  $BCF query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$BAND/chr${N}.band.vcf.gz" 2>/dev/null \
    | awk '{print $1":"$2":"$3":"$4}' | sort -u > "$OUT/.b37.$N" || return 4
  $BCF query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' "$KEYED/chr${N}.keyed38.vcf.gz" 2>/dev/null \
    | awk -v F="$OUT/.b37.$N" 'BEGIN{while((getline l < F)>0) b[l]=1}
        { k=$3; sub(/^chr/,"",k); if (k in b) print $1"\t"($2-1)"\t"$2"\t"$3"\t"$4"\t"$5 }' \
    | sort -k1,1 -k2,2n > "$OUT/.b38.$N.bed" || return 5
  local nb=$(wc -l < "$OUT/.b38.$N.bed")
  [ "$nb" -gt 0 ] || { echo "GATE FAIL chr$N: 0 band variants mapped to 38"; return 6; }
  $BED intersect -a "$OUT/.b38.$N.bed" -b "$OUT/cds.bed.gz" -u > "$OUT/.cds.$N.bed" 2>/dev/null
  local nc=$(wc -l < "$OUT/.cds.$N.bed")
  printf 'chr%s\t%s\t%s\n' "$N" "$nb" "$nc" > "$OUT/chr${N}.a.tsv"
  rm -f "$OUT/.b37.$N"
  touch "$OUT/chr${N}.a.done"
  echo "[chr$N] band38=$nb inCDS=$nc"
}
export -f one; export OUT BCF BED KEYED BAND
seq 1 22 | xargs -P ${XP:-6} -I{} bash -c 'one {}'

cat $OUT/chr*.a.tsv 2>/dev/null | sort -V > "$OUT/stageA.tsv"
NDONE=$(ls $OUT/chr*.a.done 2>/dev/null | wc -l)
awk -F'\t' -v nd="$NDONE" '{b+=$2; c+=$3} END{
  printf "STAGE_A chroms=%d band38=%d inCDS=%d (%.3f%%)\n", nd, b, c, c/b*100 }' "$OUT/stageA.tsv"
[ "$NDONE" -eq 22 ] || { echo "INCOMPLETE $NDONE/22"; exit 7; }

{ echo '##fileformat=VCFv4.2'; printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n';
  for N in $(seq 1 22); do
    awk -F'\t' '{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t.\t.\t."}' "$OUT/.cds.$N.bed" 2>/dev/null
  done | sort -k1,1V -k2,2n -u; } | gzip -c > "$OUT/vep_input.vcf.gz"
echo "VEP input: $(zcat $OUT/vep_input.vcf.gz | grep -vc '^#') variants"
echo STAGE_A_COMPLETE
