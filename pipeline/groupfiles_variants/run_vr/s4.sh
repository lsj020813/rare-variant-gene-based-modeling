: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
P2=${PLINK2:-plink2}
W=${PROJECT_ROOT}/work/run_vr; PR=${PRUNED_PLINK_PREFIX:?Set PRUNED_PLINK_PREFIX}
mkdir -p $W/gf
for N in 21 22; do
  $P2 --vcf $W/g/chr$N.vcf.gz --make-bed --const-fid 0 \
      --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 60 missing \
      --out $W/gf/chr$N > $W/gf/chr$N.log 2>&1
done
echo "rare fam: $(head -1 $W/gf/chr21.fam | awk '{print $1\" | \"$2}')"
echo "pruned fam: $(head -1 $PR.fam | awk '{print $1\" | \"$2}')"
echo $W/gf/chr22 > $W/gf/ml.txt
$P2 --bfile $W/gf/chr21 --pmerge-list $W/gf/ml.txt bfile --make-bed --out $W/gf/rare > $W/gf/r.log 2>&1
echo "rare merged: $(wc -l < $W/gf/rare.bim 2>/dev/null) markers $(wc -l < $W/gf/rare.fam 2>/dev/null) samples"
echo $W/gf/rare > $W/gf/ml2.txt
$P2 --bfile $PR --pmerge-list $W/gf/ml2.txt bfile --make-bed --out $W/gf/final > $W/gf/f.log 2>&1
if [ -s $W/gf/final.bed ]; then echo "FINAL OK: $(wc -l < $W/gf/final.bim) markers $(wc -l < $W/gf/final.fam) samples";
else echo FINAL_FAIL; grep -m3 -i 'error' $W/gf/f.log; fi
echo S4_DONE
