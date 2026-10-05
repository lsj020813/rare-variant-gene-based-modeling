#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
cd ${PROJECT_ROOT}/work/run_annot_bbj
export TMPDIR=${PROJECT_ROOT}/work/tmp TMP=${PROJECT_ROOT}/work/tmp TEMP=${PROJECT_ROOT}/work/tmp
rm -f ${PROJECT_ROOT}/work/ref/annot/extract_bbj/chr20.annot.tsv.tmp ${PROJECT_ROOT}/work/ref/annot/extract_bbj/chr20.cadd.tsv.tmp
( bash pool_bbj.sh "20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1" 6 > pool_full.log 2>&1; echo "pool_full rc=$?" >> stage4.log ) &
( bash gt_pool_bbj.sh "19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1" 2 > gt_full.log 2>&1; echo "gt_full rc=$?" >> stage4.log ) &
sleep 5
bash ${PROJECT_ROOT}/work/run_band15/wd5.sh 'ext_bb[j].sh|gpn_bb[j].sh|t2_bb[j].py|pool_bb[j].sh' ${PROJECT_ROOT}/work/run_annot_bbj 38 40 >> wd_stage4.log 2>&1 &
wait
echo STAGE4_DONE >> stage4.log
