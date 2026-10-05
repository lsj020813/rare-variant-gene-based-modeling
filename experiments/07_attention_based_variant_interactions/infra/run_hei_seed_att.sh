#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT to the working data root}"
[[ "$PROJECT_ROOT" =~ [^[:space:]] ]] || { printf "%s\n" "PROJECT_ROOT must be nonblank" >&2; exit 2; }
CODE_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
export CODE_ROOT
cd ${PROJECT_ROOT}/work/prs; PY="${PYTHON_BIN:-python3}"
startwd() { tmux kill-session -t wd_$2 2>/dev/null; tmux new-session -d -s wd_$2 "bash \"${CODE_ROOT}/infra/wd_prs_big.sh\" '$1' $2"; sleep 5
  pgrep -f "wd_prs_big.sh .*$2" >/dev/null && [ -s logs/wd_$2.log ] || { echo "WATCHDOG FAILED $2"; exit 9; }; echo "watchdog $2 ok $(date -Is)"; }
ramgate() { while [ $(free -g | awk '/Mem:/{print $7}') -lt 120 ]; do echo "waiting RAM $(date -Is)"; sleep 60; done; }
run() { ramgate; startwd '[t]rain_att_hei_seed.py' attseed_$1_$2; $PY "${CODE_ROOT}/models/attention/train_att_hei_seed.py" $1 $2 > logs/attseed_$1_$2.log 2>&1; rc=$?; tmux kill-session -t wd_attseed_$1_$2; echo "$1 $2 rc=$rc $(date -Is)"; return $rc; }
while [ ! -s out/hei_seed/extract.done ] || [ ! -s out/hei_seed/feat/feat.done ]; do sleep 60; done
run P6 smoke || { echo "P6 smoke FAILED"; exit 1; }
run P5 smoke || { echo "P5 smoke FAILED"; exit 1; }
while [ ! -s out/hei_seed/gt/extract.done ]; do sleep 60; done
run P5h smoke || { echo "P5h smoke FAILED"; exit 1; }
run P5h train || exit 2
run P6 train || exit 2
run P5 train || exit 2
echo "ok $(date -Is)" > out/hei_seed/att.done
