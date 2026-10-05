#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
OUT=${PROJECT_ROOT}/work/ref/deductive
mkdir -p "$OUT"
exec 9>"$OUT/.g5.lock"; flock -n 9 || { echo "another g5 run"; exit 0; }
BCF=${BCFTOOLS:-bcftools}
GTF=${PROJECT_ROOT}/work/ref/annotation_ref/annotation_data/gencode/gencode.v44.annotation.gtf.gz
KEYED=${PROJECT_ROOT}/work/ref/lift38_keyed
BAND=${PROJECT_ROOT}/work/ref/band_vcf

if [ ! -s "$OUT/cds.bed.gz" ]; then
  zcat "$GTF" | awk -F'\t' '$3=="CDS" && $9 ~ /gene_type "protein_coding"/ {
      match($9, /gene_id "[^"]+"/); g=substr($9, RSTART+9, RLENGTH-10);
      print $1"\t"($4-1)"\t"$5"\t"g }' | sort -k1,1 -k2,2n | gzip -c > "$OUT/cds.bed.gz"
fi
echo "CDS intervals: $(zcat $OUT/cds.bed.gz | wc -l)"

: > "$OUT/stageA.tsv"
for N in $(seq 1 22); do
  K=$KEYED/chr${N}.keyed38.vcf.gz
  [ -s "$K" ] || { echo "MISSING $K"; exit 4; }
  $BCF query -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' "$BAND/chr${N}.band.vcf.gz" 2>/dev/null \
    | awk '{print $1":"$2":"$3":"$4}' | sort -u > "$OUT/.band37.chr${N}"
  $BCF query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' "$K" 2>/dev/null \
    | awk -v F="$OUT/.band37.chr${N}" 'BEGIN{while((getline l < F)>0) b[l]=1}
        { k=$3; sub(/^chr/,"",k); if (k in b) print $1"\t"($2-1)"\t"$2"\t"$3"\t"$4"\t"$5 }' \
    | sort -k1,1 -k2,2n > "$OUT/.band38.chr${N}.bed"
  nb=$(wc -l < "$OUT/.band38.chr${N}.bed")
  ncds=$(bedtools intersect -a "$OUT/.band38.chr${N}.bed" -b "$OUT/cds.bed.gz" -u 2>/dev/null | wc -l)
  printf 'chr%s\t%s\t%s\n' "$N" "$nb" "$ncds" >> "$OUT/stageA.tsv"
  echo "[chr$N] band-on-38 $nb | in CDS $ncds"
done
awk -F'\t' '{b+=$2; c+=$3} END{printf "STAGE_A_TOTAL band=%d in_CDS=%d (%.3f%%)\n", b, c, c/b*100}' "$OUT/stageA.tsv"

{ echo '##fileformat=VCFv4.2'; echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO';
  for N in $(seq 1 22); do
    bedtools intersect -a "$OUT/.band38.chr${N}.bed" -b "$OUT/cds.bed.gz" -u 2>/dev/null \
      | awk -F'\t' '{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t.\t.\t."}'
  done; } | gzip -c > "$OUT/vep_input.vcf.gz"
echo "VEP input variants: $(zcat $OUT/vep_input.vcf.gz | grep -vc '^#')"
echo STAGE_A_COMPLETE
