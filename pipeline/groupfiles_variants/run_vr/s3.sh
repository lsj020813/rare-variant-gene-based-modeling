: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
P2=${PLINK2:-plink2}
W=${PROJECT_ROOT}/work/run_vr; mkdir -p $W/gt
for N in 21 22; do
  $P2 --vcf $W/g/chr$N.vcf.gz --make-bed --double-id \
      --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 60 missing \
      --out $W/gt/chr$N > $W/gt/chr$N.log 2>&1
  echo "chr$N markers $(wc -l < $W/gt/chr$N.bim)  first ID: $(head -1 $W/gt/chr$N.bim | cut -f2)"
done
echo $W/gt/chr22 > $W/gt/ml.txt
$P2 --bfile $W/gt/chr21 --pmerge-list $W/gt/ml.txt bfile --make-bed --out $W/gt/m > $W/gt/m.log 2>&1
if [ -s $W/gt/m.bed ]; then echo "MERGE OK: $(wc -l < $W/gt/m.bim) markers $(wc -l < $W/gt/m.fam) samples";
else echo MERGE_FAIL; grep -m2 -i error $W/gt/m.log; fi
echo S3_DONE
