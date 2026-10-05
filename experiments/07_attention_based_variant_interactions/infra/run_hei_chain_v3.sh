#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT to the working data root}"
[[ "$PROJECT_ROOT" =~ [^[:space:]] ]] || { printf "%s\n" "PROJECT_ROOT must be nonblank" >&2; exit 2; }
CODE_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
export CODE_ROOT
cd ${PROJECT_ROOT}/work/prs; PY="${PYTHON_BIN:-python3}"

startwd() {
  tmux kill-session -t wd_$3 2>/dev/null
  tmux new-session -d -s wd_$3 "bash \"${CODE_ROOT}/infra/$1\" '$2' $3"; sleep 5
  pgrep -f "$1 .*$3" >/dev/null && [ -s logs/wd_$3.log ] || { echo "WATCHDOG FAILED $3 -> abort"; exit 9; }
  echo "watchdog $3 ok $(date -Is)"
}

while [ ! -s out/hei/extract.done ]; do sleep 60; done
startwd wd_prs_big.sh '[h]ei_linear_v2.py' heilin
OMP_NUM_THREADS=32 OPENBLAS_NUM_THREADS=32 MKL_NUM_THREADS=32 $PY "${CODE_ROOT}/models/linear/hei_linear_v2.py" > logs/hei_linear.log 2>&1 || exit 3
echo "chain ok $(date -Is)" > out/hei/chain_linear.done; tmux kill-session -t wd_heilin
