: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
BCF=${BCFTOOLS:-bcftools}
P2=${PLINK2:-plink2}
W=${PROJECT_ROOT}/work/run_vr; mkdir -p $W/gs
V=${GENOTYPE_DIR:?Set GENOTYPE_DIR}/${GENOTYPE_FILENAME_PREFIX:?Set GENOTYPE_FILENAME_PREFIX}22${GENOTYPE_FILENAME_SUFFIX:?Set GENOTYPE_FILENAME_SUFFIX}
O=$W/gs/chr22
awk '{print $1"\t"$2}' $W/hc/chr22.tsv > $O.reg
$BCF view -R $O.reg -Oz -o $O.vcf.gz $V 2>$O.err
echo "records $($BCF view -H $O.vcf.gz | wc -l)"
$P2 --vcf $O.vcf.gz --make-bed --double-id --out $O > $O.p2.log 2>&1
if [ -s $O.bed ]; then
  $P2 --bfile $O --freq counts --out $O.f >/dev/null 2>&1
  awk 'NR>1{c=$6+0; if(c>10&&c<=20.5)a++; else if(c>20.5)b++; else z++} END{printf "  MAC<=10:%d  10-20.5:%d  >20.5:%d\n",z,a,b}' $O.f.acount
  echo "markers $(wc -l < $O.bim) samples $(wc -l < $O.fam)"
else echo PLINK_FAIL; tail -4 $O.p2.log; fi
echo S2_DONE
