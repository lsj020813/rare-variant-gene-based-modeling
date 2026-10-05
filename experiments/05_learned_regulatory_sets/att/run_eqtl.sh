#!/bin/bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -u
ATT=${PROJECT_ROOT}/work/fset/att
EQ=${PROJECT_ROOT}/work/ref/eqtl_cat
PY=python3
export TMPDIR=${PROJECT_ROOT}/work/tmp
cd $ATT
echo "EQTL_START $(date -Is)" | tee -a $EQ/run.log

filt () {
  nice -n 19 ionice -c3 zcat $EQ/$1.all.tsv.gz \
    | awk -F'\t' 'NR==FNR{g[$1]=1; next} FNR>1 && ($17 in g)' \
          $ATT/domain_genes.txt - \
    > $EQ/$2.filtered.tsv
  echo "$2 filtered rows=$(wc -l < $EQ/$2.filtered.tsv) $(date -Is)" \
    | tee -a $EQ/run.log
}

filt QTD000266 liver &
filt QTD000356 blood &
wait
echo "FILTER_DONE $(date -Is)" | tee -a $EQ/run.log

for t in liver blood; do
  nice -n 19 $PY eqtl_coherence.py $EQ/$t.filtered.tsv $t \
    > $EQ/$t.analysis.log 2>&1
  echo "$t analysis exit=$? $(date -Is)" | tee -a $EQ/run.log
done
touch $EQ/eqtl.done
echo "EQTL_DONE $(date -Is)" | tee -a $EQ/run.log
