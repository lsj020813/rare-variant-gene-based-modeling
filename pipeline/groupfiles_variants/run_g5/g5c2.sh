#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
O=${PROJECT_ROOT}/work/ref/deductive
L=${PROJECT_ROOT}/work/ref/loftee
V=${PROJECT_ROOT}/work/ref/vep_env
S=$V/share/ensembl-vep-116.1-0
BCF=${BCFTOOLS:-bcftools}
export PERL5LIB="$S/modules:$S:$V/lib/perl5/site_perl:${PERL5LIB:-}"
export PATH="$V/bin:$PATH"
exec 9>"$O/.g5c2.lock"; flock -n 9 || { echo "another stage C2"; exit 0; }

NAMES=$($V/bin/perl -e 'use Bio::DB::BigWig; my $w=Bio::DB::BigWig->new(-bigwig=>"'$L'/gerp.bw");
  my @s=$w->seq_ids; print join(",", @s[0..2]);' 2>/dev/null)
echo "gerp seq_ids sample: $NAMES"
case "$NAMES" in chr*) PFX="chr";; *) PFX="";; esac
echo "gerp naming: '${PFX}N' -> input must match"

if [ ! -s "$O/gencode.nochr.gtf.gz.tbi" ]; then
  zcat "$O/gencode.sorted.gtf.gz" | sed 's/^chr//' \
    | ${BGZIP:-bgzip} -@4 -c > "$O/gencode.nochr.gtf.gz"
  ${TABIX:-tabix} -p gff "$O/gencode.nochr.gtf.gz"
fi
FA=$O/GRCh38_nochr.fa
if [ ! -s "$FA.fai" ]; then
  sed 's/^>chr/>/' ${PROJECT_ROOT}/work/tmp/ref/GRCh38_no_alt.fa > "$FA"
  ${SAMTOOLS:-samtools} faidx "$FA"
fi
zcat "$O/vep_input.vcf.gz" | sed 's/^chr//' | gzip -c > "$O/vep_input_nochr.vcf.gz"
NIN=$(zcat "$O/vep_input_nochr.vcf.gz" | grep -vc '^#')
[ "$NIN" -eq 69750 ] || { echo "GATE FAIL input $NIN != 69750"; exit 4; }

$V/bin/perl $V/bin/vep -i "$O/vep_input_nochr.vcf.gz" -o "$O/vep_lof2.tsv" \
  --gtf "$O/gencode.nochr.gtf.gz" --fasta "$FA" --format vcf --tab --no_stats --force_overwrite \
  --fields "Uploaded_variation,Location,Allele,Gene,Feature,Consequence,BIOTYPE,CANONICAL,LoF,LoF_filter,LoF_flags,LoF_info" \
  --canonical --biotype --fork 6 --buffer_size 5000 --dir_plugins "$S" \
  --plugin LoF,loftee_path:$S,human_ancestor_fa:$L/human_ancestor.fa.gz,conservation_file:$L/loftee.sql,gerp_bigwig:$L/gerp.bw \
  2>&1 | grep -vE "Ignoring '" | tail -10

NL=$(grep -vc '^#' "$O/vep_lof2.tsv" 2>/dev/null || echo 0)
[ "$NL" -ge 1000 ] || { echo "GATE FAIL vep_lof2 only $NL rows"; exit 5; }

NZ=$(grep -o 'GERP_DIST:[-0-9.eE+]*' "$O/vep_lof2.tsv" | grep -vc 'GERP_DIST:0$' || echo 0)
echo "GERP_DIST non-zero rows: $NZ"
[ "$NZ" -gt 0 ] || { echo "GATE FAIL GERP still all-zero — naming/reader still broken"; exit 6; }

awk -F'\t' 'BEGIN{OFS="\t"} /^#/{next}
  $7=="protein_coding" && $8=="YES" &&
  $6 ~ /stop_gained|frameshift_variant|splice_acceptor_variant|splice_donor_variant|start_lost|transcript_ablation/ {
    print $1, $2, $4, $6, $9, $10, $12 }' "$O/vep_lof2.tsv" | sort -u > "$O/dz_all2.tsv"
awk -F'\t' '$5=="HC"' "$O/dz_all2.tsv" > "$O/dz_hc2.tsv"

echo "== candidates: $(wc -l < $O/dz_all2.tsv)"
echo "== HC:         $(wc -l < $O/dz_hc2.tsv)"
echo "== verdicts =="
awk -F'\t' '{v=($5==""||$5=="-")?"(none)":$5; c[v]++} END{for(k in c) printf "   %-8s %s\n", k, c[k]}' "$O/dz_all2.tsv" | sort -k2 -rn
echo "== LC reasons =="
awk -F'\t' '$5=="LC"{n=split($6,a,","); for(i=1;i<=n;i++) c[a[i]]++} END{for(k in c) printf "   %-18s %s\n", k, c[k]}' "$O/dz_all2.tsv" | sort -k2 -rn
echo "== 50bp rule =="
grep -o '50_BP_RULE:[A-Z]*' "$O/dz_all2.tsv" | sort | uniq -c | awk '{printf "   %s %s\n", $2, $1}'
echo "== genes =="
awk -F'\t' '{g[$3]++} END{n1=0;n2=0; for(k in g){n1++; if(g[k]>=2)n2++}
  printf "   HC genes >=1: %d\n   HC genes >=2: %d\n", n1, n2}' "$O/dz_hc2.tsv"
echo STAGE_C2_COMPLETE
