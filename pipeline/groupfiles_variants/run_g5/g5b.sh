#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
O=${PROJECT_ROOT}/work/ref/deductive
VENV=${PROJECT_ROOT}/work/ref/vep_env
VEP="$VENV/bin/perl $VENV/bin/vep"
export PERL5LIB="$VENV/share/ensembl-vep-116.1-0/modules:$VENV/share/ensembl-vep-116.1-0:$VENV/lib/perl5/site_perl:${PERL5LIB:-}"
export PATH="$VENV/bin:$PATH"
GTF=${PROJECT_ROOT}/work/ref/deductive/gencode.sorted.gtf.gz
FA=${PROJECT_ROOT}/work/tmp/ref/GRCh38_no_alt.fa
exec 9>"$O/.g5b.lock"; flock -n 9 || { echo "another stage B"; exit 0; }

test -s "$GTF.tbi" || { echo "FAIL missing tabix index for $GTF"; exit 4; }
test -s "$FA" || { echo "FAIL missing FASTA $FA"; exit 4; }

$VEP --input_file "$O/vep_input.vcf.gz" --output_file "$O/vep_out.tsv" \
  --gtf "$GTF" --fasta "$FA" --format vcf --tab --no_stats --force_overwrite \
  --fields "Uploaded_variation,Location,Allele,Gene,Feature,Consequence,BIOTYPE,CANONICAL" \
  --canonical --biotype --fork 6 --buffer_size 5000 --no_check_variants_order

NL=$(grep -vc '^#' "$O/vep_out.tsv" 2>/dev/null || echo 0)
[ "$NL" -ge 1000 ] || { echo "GATE FAIL vep_out only $NL rows"; exit 5; }
echo "VEP annotated rows: $NL"

awk -F'\t' 'BEGIN{OFS="\t"}
  /^#/ {next}
  $7=="protein_coding" && $8=="YES" {
    c=$6
    d = (c ~ /stop_gained|frameshift_variant|splice_acceptor_variant|splice_donor_variant|start_lost|transcript_ablation/) ? 1 : 0
    if (d) print $1, $2, $4, $6
  }' "$O/vep_out.tsv" | sort -u > "$O/deductive_zone.tsv"
echo "DEDUCTIVE variants: $(wc -l < $O/deductive_zone.tsv)"
awk -F'\t' '{g[$3]++} END{n1=0;n2=0; for(k in g){n1++; if(g[k]>=2) n2++}
  print "genes with >=1: " n1 "\ngenes with >=2: " n2}' "$O/deductive_zone.tsv"
awk -F'\t' '{n=split($4,a,","); for(i=1;i<=n;i++) c[a[i]]++} END{for(k in c) printf "  %-34s %s\n", k, c[k]}' \
  "$O/deductive_zone.tsv" | sort -k2 -rn
echo STAGE_B_COMPLETE
