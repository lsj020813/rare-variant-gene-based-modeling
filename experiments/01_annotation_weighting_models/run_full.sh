#!/bin/bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -o pipefail
export TMPDIR=${PROJECT_ROOT}/work/tmp TMP=${PROJECT_ROOT}/work/tmp TEMP=${PROJECT_ROOT}/work/tmp HF_HOME=${PROJECT_ROOT}/work/run_trackB/hf_cache
cd ${PROJECT_ROOT}/work/run_trackB
date >> full.log
nice -n 10 ./venv/bin/python enformer_score_v2.py --variants variants_ordered.tsv --out scores_full.tsv >> full.log 2>&1
echo EXIT $? >> full.log; date >> full.log
touch run_full.done
