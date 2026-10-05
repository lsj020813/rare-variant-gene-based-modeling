#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT to the working data root}"
[[ "$PROJECT_ROOT" =~ [^[:space:]] ]] || { printf "%s\n" "PROJECT_ROOT must be nonblank" >&2; exit 2; }
CODE_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
export CODE_ROOT
cd ${PROJECT_ROOT}/work/prs; PY="${PYTHON_BIN:-python3}"
while [ ! -s out/hei_seed/extract.done ]; do sleep 60; done
while [ $(free -g | awk '/Mem:/{print $7}') -lt 140 ]; do echo "waiting RAM $(date -Is)"; sleep 60; done
tmux kill-session -t wd_seedlin 2>/dev/null; tmux new-session -d -s wd_seedlin "bash \"${CODE_ROOT}/infra/wd_prs_big.sh\" '[h]ei_seed_linear.py' seedlin"; sleep 5
pgrep -f 'wd_prs_big.sh .*seedlin' >/dev/null && [ -s logs/wd_seedlin.log ] || { echo "WATCHDOG FAILED seedlin"; exit 9; }; echo "watchdog seedlin ok $(date -Is)"
OMP_NUM_THREADS=32 OPENBLAS_NUM_THREADS=32 MKL_NUM_THREADS=32 $PY "${CODE_ROOT}/models/linear/hei_seed_linear.py" > logs/hei_seed_linear.log 2>&1; rc=$?
tmux kill-session -t wd_seedlin; echo "linear rc=$rc $(date -Is)"; [ $rc -eq 0 ] && echo "ok $(date -Is)" > out/hei_seed/linear.done
