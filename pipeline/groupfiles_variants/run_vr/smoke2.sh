: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
BCF=${BCFTOOLS:-bcftools}
P2=${PLINK2:-plink2}
W=${PROJECT_ROOT}/work/run_vr; mkdir -p $W/geno_smoke
V=${GENOTYPE_DIR:?Set GENOTYPE_DIR}/${GENOTYPE_FILENAME_PREFIX:?Set GENOTYPE_FILENAME_PREFIX}22${GENOTYPE_FILENAME_SUFFIX:?Set GENOTYPE_FILENAME_SUFFIX}
O=$W/geno_smoke/chr22
awk '{print $1"\t"$2}' $W/chunks/chr22.rare.tsv > $O.regions
S=$(date +%s)
$BCF view -R $O.regions -Oz -o $O.vcf.gz $V 2>$O.err
echo "extract $(( $(date +%s) - S ))s  records $($BCF view -H $O.vcf.gz | wc -l)"
$P2 --vcf $O.vcf.gz dosage=DS --make-bed --double-id --out $O > $O.p2.log 2>&1
if [ -s $O.bed ]; then echo "plink markers $(wc -l < $O.bim) samples $(wc -l < $O.fam)";
  $P2 --bfile $O --freq counts --out $O.f >/dev/null 2>&1
  awk 'NR>1{c=$6+0; if(c>10&&c<=20.5)a++; else if(c>20.5)b++; else z++} END{printf "  MAC<=10:%d  10-20.5:%d  >20.5:%d\n",z,a,b}' $O.f.acount
else echo 'PLINK FAILED'; tail -5 $O.p2.log; fi
echo SMOKE2_DONE
