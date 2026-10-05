#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT to the working data root}"
[[ "$PROJECT_ROOT" =~ [^[:space:]] ]] || { printf "%s\n" "PROJECT_ROOT must be nonblank" >&2; exit 2; }
CODE_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
export CODE_ROOT
cd ${PROJECT_ROOT}/work/prs; PY="${PYTHON_BIN:-python3}"
startwd() { tmux kill-session -t wd_$2 2>/dev/null; tmux new-session -d -s wd_$2 "bash \"${CODE_ROOT}/infra/wd_prs_big.sh\" '$1' $2"; sleep 5
  pgrep -f "wd_prs_big.sh .*$2" >/dev/null && [ -s logs/wd_$2.log ] || { echo "WATCHDOG FAILED $2"; exit 9; }; echo "watchdog $2 ok $(date -Is)"; }
ramgate() { while [ $(free -g | awk '/Mem:/{print $7}') -lt 120 ]; do echo "waiting RAM $(date -Is)"; sleep 60; done; }
run() { ramgate; startwd '[t]rain_att_hei_res_v4[.]py' attres_$1_$2; $PY "${CODE_ROOT}/models/attention/train_att_hei_res_v4.py" $1 $2 > logs/attres_$1_$2.log 2>&1; rc=$?; tmux kill-session -t wd_attres_$1_$2; echo "$1 $2 rc=$rc $(date -Is)"; return $rc; }
while [ ! -s out/hei_seed/p3ann128_oof.json ]; do sleep 60; done
echo "p3 oof ready $(date -Is)"
ramgate; startwd '[t]rain_att_hei_res_v4[.]py' attres_logtest; HEI_LOGTEST=1 $PY "${CODE_ROOT}/models/attention/train_att_hei_res_v4.py" P5hres train > logs/attres_P5hres_logtest.log 2>&1; rc=$?; tmux kill-session -t wd_attres_logtest; echo "logtest rc=$rc $(date -Is)"; [ $rc -eq 0 ] || { echo "LOGTEST FAILED"; exit 1; }
run P6res smoke || { echo "P6res smoke FAILED"; exit 1; }
run P5res smoke || { echo "P5res smoke FAILED"; exit 1; }
run P5hres smoke || { echo "P5hres smoke FAILED"; exit 1; }
run P6res train || exit 2
run P5hres train || exit 2
run P5res train || exit 2
echo "ok $(date -Is)" > out/hei_seed/att_res.done
