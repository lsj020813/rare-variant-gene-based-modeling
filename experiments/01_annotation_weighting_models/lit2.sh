: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
D=${PROJECT_ROOT}/work/ref/saige_step2_bwg
cat $D/{dm,tchl,htn,lip}.chr[0-9]* 2>/dev/null | awk -F'\t' 'FNR==1{next} {print FILENAME"\t"$0}' > /dev/null
for t in dm tchl htn lip; do for f in $D/$t.chr[0-9]*; do awk -F'\t' -v t=$t 'NR>1{print t"\t"$1"\t"$4"\t"$12}' $f; done; done > all_truth_p.tsv
while read id s; do printf '%-8s ' $s; for t in dm tchl htn lip; do p=$(awk -F'\t' -v t=$t -v g="$id" '$1==t && $2==g{printf "%.2e(m=%s)", $3, $4; exit}' all_truth_p.tsv); printf '%s=%s ' $t "${p:-NA}"; done; echo; done < lit_genes.map
echo LIT2_DONE
