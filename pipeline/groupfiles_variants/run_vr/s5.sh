: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
P2=${PLINK2:-plink2}
W=${PROJECT_ROOT}/work/run_vr; PR=${PRUNED_PLINK_PREFIX:?Set PRUNED_PLINK_PREFIX}
$P2 --bfile $PR            --make-pgen --sort-vars --out $W/gf/pr_s  > $W/gf/a.log 2>&1
$P2 --bfile $W/gf/rare     --make-pgen --sort-vars --out $W/gf/ra_s  > $W/gf/b.log 2>&1
echo "pruned pgen: $([ -s $W/gf/pr_s.pgen ] && echo Y || echo N)  rare pgen: $([ -s $W/gf/ra_s.pgen ] && echo Y || echo N)"
echo $W/gf/ra_s > $W/gf/ml3.txt
$P2 --pfile $W/gf/pr_s --pmerge-list $W/gf/ml3.txt pfile --make-bed --out $W/gf/final2 > $W/gf/c.log 2>&1
if [ -s $W/gf/final2.bed ]; then
  echo "MERGE OK: $(wc -l < $W/gf/final2.bim) markers $(wc -l < $W/gf/final2.fam) samples"
  $P2 --bfile $W/gf/final2 --freq counts --out $W/gf/final2f >/dev/null 2>&1
  awk 'NR>1{c=$6+0; if(c>10&&c<=20.5)a++; else if(c>20.5)b++; else z++} END{printf "  MAC<=10:%d  10-20.5:%d  >20.5:%d\n",z,a,b}' $W/gf/final2f.acount
else echo MERGE_FAIL; grep -m3 -i error $W/gf/c.log; fi
echo S5_DONE
