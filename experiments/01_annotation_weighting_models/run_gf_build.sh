#!/usr/bin/env bash
: "${CONDA_INIT_SCRIPT:?Set CONDA_INIT_SCRIPT}"
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -o pipefail
W=${PROJECT_ROOT}/work; B=$W/run_l3b
export TMPDIR=$W/tmp TMP=$W/tmp TEMP=$W/tmp; mkdir -p "$TMPDIR" "$B/out"
source ${CONDA_INIT_SCRIPT} && conda activate "${CONDA_ENV:?Set CONDA_ENV}"
cd "$B"
nice -n 19 ionice -c3 python3 mk_gf_chr19.py \
  --root $W/ref --cadd-dir $W/ref/features --cadd-source features --chrom 19 \
  --group-file $W/ref/groupfiles_bwg/chr19.B_3kb_re2g.txt \
  --phi-json $W/run_band15/model_v10_out/e6_runs/primary_f0/phi.json \
  --v10-script $W/run_band15/model_v10_out/l1_train_v10.py \
  --t1-partial-policy blank --missing-phi-rule intercept-prior \
  --threads 2 --out $B/out
echo "GF_BUILD_RC=$?"
