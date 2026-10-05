#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT to the working data root}"
[[ "$PROJECT_ROOT" =~ [^[:space:]] ]] || { printf "%s\n" "PROJECT_ROOT must be nonblank" >&2; exit 2; }
set -uo pipefail
N=$1; W=${PROJECT_ROOT}/work/prs; B="${BCFTOOLS_BIN:-bcftools}"; P="${PLINK2_BIN:-plink2}"
mkdir -p $W/geno_hm3; O=$W/geno_hm3/chr$N
[ -s $O.done ] && { echo "have chr$N"; exit 0; }
awk -v c=$N 'NR>1 && $1==c {print $1"\t"$3}' $W/ref/ldblk_1kg_eas/snpinfo_1kg_hm3 > $O.pos
$B view -R $O.pos -Ou ${PROJECT_ROOT}/work/ref/orig_index/chr$N.vcf.gz 2>/dev/null \
 | $P --bcf /dev/stdin dosage=DS --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 200 --make-pgen --threads 2 --memory 8000 --out $O.tmp > $O.plink.log 2>&1
rc=$?
if [ $rc -eq 0 ] && [ -s $O.tmp.pgen ]; then for x in pgen pvar psam; do mv $O.tmp.$x $O.$x; done; echo "ok $(date -Is) $(wc -l < $O.pvar)" > $O.done; fi
echo "chr$N rc=$rc"
