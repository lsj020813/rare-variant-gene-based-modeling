: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
P2=${PLINK2:-plink2}
W=${PROJECT_ROOT}/work/run_vr; PR=${PRUNED_PLINK_PREFIX:?Set PRUNED_PLINK_PREFIX}
echo $W/gf/rare > $W/gf/ml4.txt
echo '--- A: --sample-inner-join ---'
$P2 --bfile $PR --pmerge-list $W/gf/ml4.txt bfile --sample-inner-join \
    --make-bed --out $W/gf/tA > $W/gf/tA.log 2>&1
[ -s $W/gf/tA.bed ] && echo "  OK $(wc -l < $W/gf/tA.bim) markers $(wc -l < $W/gf/tA.fam) samples" \
  || { echo '  FAIL'; grep -m2 -i error $W/gf/tA.log; }
echo '--- B: --sample-inner-join --pheno-inner-join ---'
$P2 --bfile $PR --pmerge-list $W/gf/ml4.txt bfile --sample-inner-join --pheno-inner-join \
    --make-bed --out $W/gf/tB > $W/gf/tB.log 2>&1
[ -s $W/gf/tB.bed ] && echo "  OK $(wc -l < $W/gf/tB.bim) markers $(wc -l < $W/gf/tB.fam) samples" \
  || { echo '  FAIL'; grep -m2 -i error $W/gf/tB.log; }
echo S7_DONE
