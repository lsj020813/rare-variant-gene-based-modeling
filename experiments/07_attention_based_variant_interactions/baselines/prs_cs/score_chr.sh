#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT to the working data root}"
[[ "$PROJECT_ROOT" =~ [^[:space:]] ]] || { printf "%s\n" "PROJECT_ROOT must be nonblank" >&2; exit 2; }
set -uo pipefail
N=$1; T=$2; W=${PROJECT_ROOT}/work/prs; P="${PLINK2_BIN:-plink2}"
mkdir -p $W/score_$T; O=$W/score_$T/chr$N
[ -s $O.done ] && exit 0
f=$(ls $W/prscs_$T/${T}_pst_eff_a1_b0.5_phiauto_chr$N.txt)
awk 'NR==FNR{m[$1]=$2;next} ($2 in m){print m[$2]"\t"$4"\t"$6}' $W/geno_hm3/chr$N.rsmap $f > $O.w
cut -f1 $W/private/split.tsv | tail -n+2 | cat > $O.keep
$P --pfile $W/geno_hm3/chr$N --keep $O.keep --score $O.w 1 2 3 cols=+scoresums --threads 2 --memory 8000 --out $O > $O.log 2>&1
rc=$?; [ $rc -eq 0 ] && [ -s $O.sscore ] && echo "ok $(wc -l < $O.w)" > $O.done; rm -f $O.keep
echo "score $T chr$N rc=$rc"
