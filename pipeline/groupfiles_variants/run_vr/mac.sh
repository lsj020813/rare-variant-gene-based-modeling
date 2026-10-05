#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
P2=${PLINK2:-plink2}
SRC=${PRUNED_PLINK_PREFIX:?Set PRUNED_PLINK_PREFIX}
W=${PROJECT_ROOT}/work/run_vr; mkdir -p $W; cd $W
$P2 --bfile "$SRC" --freq counts --out prune110k >/dev/null 2>&1
awk 'NR>1{c=$6+0; if(c>10 && c<=20.5) a++; else if(c>20.5) b++; else z++}
     END{printf "MAC<=10: %d\nMAC 10-20.5: %d\nMAC>20.5: %d\ntotal: %d\n", z,a,b,z+a+b}'     prune110k.acount
echo "MAC_DONE"
