: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
P1=${PROJECT_ROOT}/work/ref/plink19_env/bin/plink
P2=${PLINK2:-plink2}
W=${PROJECT_ROOT}/work/run_vr
PR=${PRUNED_PLINK_PREFIX:?Set PRUNED_PLINK_PREFIX}
OUT=${PROJECT_ROOT}/work/ref/vr_plink_v3; mkdir -p $OUT
: > $W/ml.txt; for N in $(seq 2 22); do echo "$W/g/chr${N}" >> $W/ml.txt; done
$P1 --bfile $W/g/chr1 --merge-list $W/ml.txt --keep-allele-order --make-bed --out $W/rare_all2 > $W/m1b.log 2>&1
[ -s $W/rare_all2.bed ] || { echo 'GATE FAIL rare merge'; grep -m3 -i error $W/m1b.log; exit 6; }
SUMIN=0; for N in $(seq 1 22); do SUMIN=$(( SUMIN + $(wc -l < $W/g/chr$N.bim) )); done
NR_=$(wc -l < $W/rare_all2.bim)
echo "rare merged: $NR_ markers (input sum $SUMIN)  samples $(wc -l < $W/rare_all2.fam)"
$P1 --bfile $PR --bmerge $W/rare_all2 --keep-allele-order --make-bed --out $OUT/vr_v3 > $W/m2b.log 2>&1
[ -s $OUT/vr_v3.bed ] || { echo 'GATE FAIL final merge'; grep -m3 -i error $W/m2b.log; exit 6; }
NM=$(wc -l < $OUT/vr_v3.bim); NS=$(wc -l < $OUT/vr_v3.fam)
EXP=$(( $(wc -l < $PR.bim) + NR_ ))
echo "FINAL: $NM markers  $NS samples   (expect $(wc -l < $PR.bim) + $NR_ = $EXP)"
[ "$NM" -eq "$EXP" ] || { echo "GATE FAIL: identity $NM != $EXP"; exit 6; }
[ "$NS" -eq ${N_SAMPLES} ] || { echo "GATE FAIL: samples $NS != ${N_SAMPLES}"; exit 6; }
DUP=$(awk '{print $2}' $OUT/vr_v3.bim | sort | uniq -d | wc -l)
[ "$DUP" -eq 0 ] || { echo "GATE FAIL: $DUP duplicate IDs"; exit 6; }
$P2 --bfile $OUT/vr_v3 --freq counts --out $W/ffz >/dev/null 2>&1
awk 'NR>1{c=$6+0; if(c>10&&c<=20.5)a++; else if(c>20.5)b++; else z++} END{printf "MAC<=10:%d  10-20.5:%d  >20.5:%d\n",z,a,b; if(a<200){print "GATE FAIL cat1"; exit 6} if(b<200){print "GATE FAIL cat2"; exit 6}}' $W/ffz.acount || exit 6
echo VR_PLINK_V3_COMPLETE
