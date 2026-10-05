#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
export TMPDIR=${PROJECT_ROOT}/work/tmp
cd ${PROJECT_ROOT}/work/run_band15/model_v8_out
exec python3 ./build_resid_v8.py --threads 4 --memory-gb 60
