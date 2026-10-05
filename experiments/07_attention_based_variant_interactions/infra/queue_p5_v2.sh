#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT to the working data root}"
[[ "$PROJECT_ROOT" =~ [^[:space:]] ]] || { printf "%s\n" "PROJECT_ROOT must be nonblank" >&2; exit 2; }
CODE_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
export CODE_ROOT
cd ${PROJECT_ROOT}/work/prs
while [ ! -s out/p6_bbj.done ]; do sleep 60; done
[ -s out/train_bbj_P6.json ] || { echo "P6 json missing" > out/p5_bbj.fail; exit 1; }
tmux new-session -d -s prs_wd10 "bash \"${CODE_ROOT}/infra/wd_prs_big.sh\" '[t]rain_att_v4.py' p5run"
sleep 2
"${PYTHON_BIN:-python3}" "${CODE_ROOT}/models/attention/train_att_v4.py" bbj P5 > logs/train_bbj_P5.log 2>&1
echo "done $(date -Is)" > out/p5_bbj.done
