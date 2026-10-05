#!/bin/bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
cd ${PROJECT_ROOT}/work/haplotype_context_20260922_v1/run_20260922T052403Z
export PYTHONPATH=${PROJECT_ROOT}/work/haplotype_context_20260922_v1/run_20260922T052403Z/code OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=python3
scan(){ for ch in $@; do
  nice -n 19 ionice -c3 $PY -c "
import json,sys
from hcf import tiles
print(json.dumps(tiles.scan_chrom('$ch','${PROJECT_ROOT}/work/ref/lift38_keyed','${PROJECT_ROOT}/work/haplotype_context_20260922_v1/run_20260922T052403Z/inventory/coding_splice_hg38.npz','${PROJECT_ROOT}/work/haplotype_context_20260922_v1/run_20260922T052403Z/variants')))
" >> ${PROJECT_ROOT}/work/haplotype_context_20260922_v1/run_20260922T052403Z/logs/inv_$1.log 2>&1
done; }
scan 1 5 9 13 17 21 &
scan 2 6 10 14 18 22 &
scan 3 7 11 15 19 &
scan 4 8 12 16 20 &
wait
echo done > ${PROJECT_ROOT}/work/haplotype_context_20260922_v1/run_20260922T052403Z/logs/inventory.done
