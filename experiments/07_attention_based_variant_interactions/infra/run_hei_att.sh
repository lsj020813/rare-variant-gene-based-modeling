#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT to the working data root}"
[[ "$PROJECT_ROOT" =~ [^[:space:]] ]] || { printf "%s\n" "PROJECT_ROOT must be nonblank" >&2; exit 2; }
CODE_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
export CODE_ROOT
cd ${PROJECT_ROOT}/work/prs; PY="${PYTHON_BIN:-python3}"
while [ ! -s private/hei/base_hei.npz ] || [ ! -s out/hei/extract.done ] || [ ! -s out/hei/feat/feat.done ]; do sleep 60; done
for A in P6 P5; do
  $PY "${CODE_ROOT}/models/attention/train_att_hei.py" $A bench > logs/bench_hei_$A.log 2>&1 || { echo "bench $A failed" ; exit 1; }
  $PY "${CODE_ROOT}/models/attention/train_att_hei.py" $A > logs/train_hei_$A.log 2>&1 || { echo "train $A failed"; exit 2; }
  echo "$A done $(date -Is)"
done
echo ok > out/hei/att.done
