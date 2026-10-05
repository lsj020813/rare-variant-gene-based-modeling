#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
export TMPDIR=${PROJECT_ROOT}/work/tmp
LOG=$1; shift
cd ${PROJECT_ROOT}/work/run_band15/model_v8_out
exec python3 ./l1_train_v8.py "$@" > ${PROJECT_ROOT}/work/run_band15/$LOG 2>&1
