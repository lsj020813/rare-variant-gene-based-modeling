#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
cd ${PROJECT_ROOT}/work/run_annot_bbj
export TMPDIR=${PROJECT_ROOT}/work/tmp
PY=python3
for i in $(seq 1 60); do n=$(ls ${PROJECT_ROOT}/work/ref/annot/t2_bbj/chr{15,16,17,18,19,20}.t2.tsv.done 2>/dev/null | wc -l); [ "$n" -eq 6 ] && break; sleep 20; done
( ulimit -v 24000000; $PY assemble_bbj.py --chrs 20,19,18,17,16,15 --out bbj_annot_sub.tsv.gz > asm_sub.log 2>&1; echo "asm_sub rc=$?" >> stage5.log ) 
[ -s ${PROJECT_ROOT}/work/ref/annot/bbj_l1/bbj_annot_sub.tsv.gz.done ] && ( ulimit -v 16000000; $PY pip_bbj_B.py --annot bbj_annot_sub.tsv.gz --outdir bbj_pip_sub --traits-out traits_sub.tsv > pipB_sub.log 2>&1; echo "pipB_sub rc=$?" >> stage5.log )
echo SUB_DONE >> stage5.log
