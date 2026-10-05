#!/bin/bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
cd ${PROJECT_ROOT}/work/haplotype_context_20260922_v1/run_20260922T052403Z
export PYTHONPATH=${PROJECT_ROOT}/work/haplotype_context_20260922_v1/run_20260922T052403Z/code OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
ulimit -v 12582912
PY=python3
date -u +'%Y-%m-%dT%H:%M:%SZ start' >> logs/build_states.log
cat /proc/loadavg >> logs/build_states.log
nice -n 19 ionice -c3 $PY -m hcf.workflow_ext build-states \
  --config config/resolved_config.yaml \
  --tilelist manifests/structure_tiles.tsv \
  --vardir variants --split private/sample_split_manifest.private.tsv \
  --outdir states --shard-dir shards --tmpdir tmp --workers 4 >> logs/build_states.log 2>&1
rc=$?
date -u +'%Y-%m-%dT%H:%M:%SZ end' >> logs/build_states.log
echo $rc > logs/build_states.done
