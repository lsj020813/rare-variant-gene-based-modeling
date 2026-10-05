: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
BCF=${BCFTOOLS:-bcftools}
V=${GENOTYPE_DIR:?Set GENOTYPE_DIR}/${GENOTYPE_FILENAME_PREFIX:?Set GENOTYPE_FILENAME_PREFIX}22${GENOTYPE_FILENAME_SUFFIX:?Set GENOTYPE_FILENAME_SUFFIX}
W=${PROJECT_ROOT}/work/run_vr; S=$(date +%s)
$BCF query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/R2\t[%GT]\n' $V 2>/dev/null | awk 'BEGIN{n=0}{ if($5+0<0.3) next; s=$6; ac=gsub(/1/,"1",s);   if(ac>10 && ac<=20){print $1"\t"$2"\t"$3"\t"$4"\t"ac; n++; if(n>=50) exit} }' > $W/hcsmoke.tsv
echo "elapsed $(( $(date +%s) - S ))s  markers $(wc -l < $W/hcsmoke.tsv)"
awk '{print $5}' $W/hcsmoke.tsv | sort -n | uniq -c | head -12
echo HCSMOKE_DONE
