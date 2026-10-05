#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT to the working data root}"
[[ "$PROJECT_ROOT" =~ [^[:space:]] ]] || { printf "%s\n" "PROJECT_ROOT must be nonblank" >&2; exit 2; }
CODE_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
export CODE_ROOT
cd ${PROJECT_ROOT}/work/prs; PY="${PYTHON_BIN:-python3}"
while [ ! -s out/hei/prscs.done ]; do sleep 60; done
seq 1 22 | xargs -P 6 -I{} bash "${CODE_ROOT}/baselines/prs_cs/score_chr.sh" {} hei
n=$(ls score_hei/*.done | wc -l); echo "scored $n $(date -Is)"; [ $n -eq 22 ] || exit 1
$PY "${CODE_ROOT}/data_prep/hei_prep_model.py" > logs/hei_prep_model.log 2>&1 || exit 2
while [ ! -s out/hei/extract.done ]; do sleep 60; done
OMP_NUM_THREADS=32 OPENBLAS_NUM_THREADS=32 MKL_NUM_THREADS=32 $PY "${CODE_ROOT}/models/linear/hei_linear_v2.py" > logs/hei_linear.log 2>&1 || exit 3
echo "chain ok $(date -Is)" > out/hei/chain_linear.done
