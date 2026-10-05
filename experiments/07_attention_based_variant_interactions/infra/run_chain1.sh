#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT to the working data root}"
[[ "$PROJECT_ROOT" =~ [^[:space:]] ]] || { printf "%s\n" "PROJECT_ROOT must be nonblank" >&2; exit 2; }
CODE_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
export CODE_ROOT
cd ${PROJECT_ROOT}/work/prs
PY="${PYTHON_BIN:-python3}"
seq 1 22 | xargs -P 11 -I{} bash "${CODE_ROOT}/baselines/prs_cs/score_chr.sh" {} bbj
n=$(ls score_bbj/*.done | wc -l); echo "scored $n"; [ $n -eq 22 ] || exit 1
$PY "${CODE_ROOT}/data_prep/prep_model.py" bbj || exit 2
$PY "${CODE_ROOT}/models/linear/train_linear.py" bbj > logs/train_linear_bbj.log 2>&1 || exit 3
echo "chain1 ok $(date -Is)" > out/chain1.done
