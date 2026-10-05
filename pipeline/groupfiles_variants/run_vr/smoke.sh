set -uo pipefail
BCF=${BCFTOOLS:-bcftools}
V=${GENOTYPE_DIR:?Set GENOTYPE_DIR}/${GENOTYPE_FILENAME_PREFIX:?Set GENOTYPE_FILENAME_PREFIX}22${GENOTYPE_FILENAME_SUFFIX:?Set GENOTYPE_FILENAME_SUFFIX}
S=$(date +%s)
$BCF query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/MAF\t%INFO/R2\n' $V 2>/dev/null | awk -v lo=${MAF_LO} -v hi=${MAF_HI} -v r2=0.3 -v want=200 'BEGIN{n=0}{m=$5+0;r=$6+0; if(m>lo&&m<=hi&&r>=r2){print $1"\t"$2"\t"$3"\t"$4; n++; if(n>=want) exit}}' > smoke_rare.tsv
echo "elapsed $(( $(date +%s) - S ))s  markers $(wc -l < smoke_rare.tsv)"
head -2 smoke_rare.tsv
echo SMOKE_DONE
