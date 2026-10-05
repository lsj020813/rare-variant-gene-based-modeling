: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
P1=${PROJECT_ROOT}/work/ref/plink19_env/bin/plink
P2=${PLINK2:-plink2}
W=${PROJECT_ROOT}/work/run_vr; PR=${PRUNED_PLINK_PREFIX:?Set PRUNED_PLINK_PREFIX}
mkdir -p $W/cmp
echo '=== plink1.9 --bmerge (rare 2chr smoke set + pruned) ==='
S=$(date +%s)
$P1 --bfile $PR --bmerge $W/gf/rare --keep-allele-order --make-bed --out $W/cmp/p19 > $W/cmp/p19.log 2>&1
if [ -s $W/cmp/p19.bed ]; then
  echo "  OK  $(( $(date +%s) - S ))s  markers $(wc -l < $W/cmp/p19.bim)  samples $(wc -l < $W/cmp/p19.fam)"
  $P2 --bfile $W/cmp/p19 --freq counts --out $W/cmp/p19f >/dev/null 2>&1
  awk 'NR>1{c=$6+0; if(c>10&&c<=20.5)a++; else if(c>20.5)b++; else z++} END{printf "  MAC<=10:%d  10-20.5:%d  >20.5:%d\n",z,a,b}' $W/cmp/p19f.acount
  echo "  dup IDs: $(awk '{print $2}' $W/cmp/p19.bim | sort | uniq -d | wc -l)"
else echo '  FAIL'; grep -m3 -i error $W/cmp/p19.log; fi
echo '=== 검산: 입력 마커 합 ==='
echo "  pruned $(wc -l < $PR.bim) + rare $(wc -l < $W/gf/rare.bim) = $(( $(wc -l < $PR.bim) + $(wc -l < $W/gf/rare.bim) ))"
echo CMP_DONE
