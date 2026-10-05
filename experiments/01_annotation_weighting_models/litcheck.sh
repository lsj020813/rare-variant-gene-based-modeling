#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
G=${PROJECT_ROOT}/work/ref/deductive/gencode.nochr.gtf.gz; D=${PROJECT_ROOT}/work/ref/saige_step2_bwg
BCF=bcftools
zcat $G | awk -F'\t' '$3=="gene"' | grep -E 'gene_name "(SLC30A8|PSMB8|TAP1|HMGCR|APOE|LDLR|CETP|PCSK9)"' | sed -E 's/.*gene_id "([^"]+)".*gene_name "([^"]+)".*/\1 \2/' > lit_genes.map
echo "--- gene p in our truth set (SKAT-O Pvalue)"
while read id s; do printf '%-8s %-22s ' $s $id; for t in dm tchl; do p=$(cat $D/$t.chr*.txt 2>/dev/null | awk -v g="$id" '$1==g{print $2; exit}'); printf '%s=%s ' $t "${p:-NA}"; done; echo; done < lit_genes.map
echo "--- SLC30A8 I349F (rs13266634 region, GRCh37 chr8:118184783) in our band VCF"
$BCF view -H -r 8:118184783 ${PROJECT_ROOT}/work/ref/band_vcf/chr8.band.vcf.gz 2>/dev/null | cut -f1-5,8 | cut -c1-200
echo "(empty above = not in 0.1-1% band)"; $BCF view -H -r 8:118184783 ${PROJECT_ROOT}/work/ref/orig_index/chr8.vcf.gz 2>/dev/null | cut -f1-5,8 | cut -c1-200 | sed 's/^/orig: /'
echo LIT_DONE
