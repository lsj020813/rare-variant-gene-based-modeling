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

while [ ! -s private/hei/base_hei.npz ] || [ ! -s out/hei/extract.done ] || [ ! -s out/hei/feat/feat.done ]; do sleep 60; done
while [ $(free -g | awk '/Mem:/{print $7}') -lt 120 ]; do echo "waiting RAM $(date -Is)"; sleep 60; done
for A in P6 P5; do
  startwd wd_prs_big.sh '[t]rain_att_hei.py' heiatt$A
  $PY "${CODE_ROOT}/models/attention/train_att_hei.py" $A bench > logs/bench_hei_$A.log 2>&1 || { echo "bench $A failed"; exit 1; }
  $PY "${CODE_ROOT}/models/attention/train_att_hei.py" $A > logs/train_hei_$A.log 2>&1 || { echo "train $A failed"; exit 2; }
  tmux kill-session -t wd_heiatt$A; echo "$A done $(date -Is)"
done
echo ok > out/hei/att.done
