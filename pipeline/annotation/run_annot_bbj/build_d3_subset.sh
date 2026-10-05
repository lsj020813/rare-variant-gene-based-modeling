#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -u
cd ${PROJECT_ROOT}/work/run_annot_bbj
PY=python3
L1=${PROJECT_ROOT}/work/ref/annot/bbj_l1
[ -s $L1/bbj_annot_d3.tsv.gz.done ] || {
  zcat $L1/bbj_annot.tsv.gz | awk -F'\t' 'NR==1 || $2==7 || $2==13 || $2==16 || $2==19' | gzip -1 > $L1/bbj_annot_d3.tsv.gz.tmp || exit 3
  N=$(( $(zcat $L1/bbj_annot_d3.tsv.gz.tmp | wc -l) - 1 )); [ "$N" -gt 0 ] || { echo GATE FAIL 0 rows; exit 4; }
  mv $L1/bbj_annot_d3.tsv.gz.tmp $L1/bbj_annot_d3.tsv.gz; echo "ok $N" > $L1/bbj_annot_d3.tsv.gz.done; echo "[d3] annot rows $N"
}
$PY pip_bbj_B.py --annot bbj_annot_d3.tsv.gz --outdir bbj_pip_d3 --traits-out traits_d3.tsv --fold-map '{"13":0,"16":1,"7":2,"19":3}' > pipB_d3.log 2>&1 || { echo GATE FAIL pipB; exit 5; }
tail -n 2 pipB_d3.log
echo D3_SUBSET_DONE
