: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
P2=${PLINK2:-plink2}
BCF=${BCFTOOLS:-bcftools}
W=${PROJECT_ROOT}/work/run_vr
PR=${PRUNED_PLINK_PREFIX:?Set PRUNED_PLINK_PREFIX}
V=${GENOTYPE_DIR:?Set GENOTYPE_DIR}/${GENOTYPE_FILENAME_PREFIX:?Set GENOTYPE_FILENAME_PREFIX}22${GENOTYPE_FILENAME_SUFFIX:?Set GENOTYPE_FILENAME_SUFFIX}
mkdir -p $W/sp
awk '$1==22{print $1"\t"$4}' $PR.bim > $W/sp/common22.reg
awk '{print $1"\t"$2}' $W/hc/chr22.tsv > $W/sp/rare22.reg
cat $W/sp/common22.reg $W/sp/rare22.reg | sort -k1,1 -k2,2n -u > $W/sp/both22.reg
echo "common $(wc -l < $W/sp/common22.reg)  rare $(wc -l < $W/sp/rare22.reg)  union $(wc -l < $W/sp/both22.reg)"
$BCF view -R $W/sp/both22.reg -Oz -o $W/sp/both22.vcf.gz $V 2>$W/sp/e.log
echo "extracted records: $($BCF view -H $W/sp/both22.vcf.gz | wc -l)"
$P2 --vcf $W/sp/both22.vcf.gz --make-bed --const-fid 0 \
    --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 60 missing \
    --out $W/sp/both22 > $W/sp/p.log 2>&1
if [ -s $W/sp/both22.bed ]; then
  echo "plink: $(wc -l < $W/sp/both22.bim) markers $(wc -l < $W/sp/both22.fam) samples"
  $P2 --bfile $W/sp/both22 --freq counts --out $W/sp/f >/dev/null 2>&1
  awk 'NR>1{c=$6+0; if(c>10&&c<=20.5)a++; else if(c>20.5)b++; else z++} END{printf "  MAC<=10:%d  10-20.5:%d  >20.5:%d\n",z,a,b}' $W/sp/f.acount
else echo PLINK_FAIL; tail -3 $W/sp/p.log; fi
echo S6_DONE
