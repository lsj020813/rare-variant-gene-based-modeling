#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
N=$1; C=${PROJECT_ROOT}/work/ref/l1/B/.claim.chr$N
mkdir $C 2>/dev/null || exit 0
python3 ${PROJECT_ROOT}/work/run_l1/bagg.py $N
