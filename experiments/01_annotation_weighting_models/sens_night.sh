#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
D=${PROJECT_ROOT}/work/run_band15/model_v10_out; R=$D/e6_runs; A=${PROJECT_ROOT}/work/ref/annot
PY=python3
L=$R/sens_night.log; say(){ echo "[sens] $(date '+%m-%d %H:%M') $*" | tee -a $L; }
export TMPDIR=${PROJECT_ROOT}/work/tmp
cd $D

common=( --fold 0 --bbj-annot $A/bbj_l1/bbj_annot.tsv.gz --bbj-pip $A/bbj_l1/bbj_pip
         --traits $A/bbj_l1/traits.tsv --annot-stats $A/cache/fm_all.stats.json
         --train-bands all --inner-va-folds 2 --lam-divisors 10 100 1000 --null-perm 50
         --inner-seed 20260908 --test-seed 123 --inner-maxiter 30 --maxiter 60 --design-chunk 16384
         --c1-columns v9-eight --alpha 0.05 --power-target 0.8 --memory-gb 60 )

CH=$($PY - <<'PYEOF'
import os as _env_os, re as _env_re
def _env_path(value):
    return _env_re.sub(r"\$\{([A-Z][A-Z0-9_]*)\}", lambda match: _env_os.environ[match.group(1)], value)
import re
s=open(_env_path('${PROJECT_ROOT}/work/run_band15/model_v10_out/l1_train_v10.py')).read()
m=re.search(r'--phi-columns[^)]{0,400}?choices\s*=\s*(\[[^\]]*\])', s, re.S)
print(m.group(1) if m else 'NOCHOICES')
PYEOF
)
say "phi-columns 선택지: $CH"

run(){ local TAG=$1; shift
  local O=$R/$TAG
  [ -s $O/L1_DONE ] && { say "$TAG 이미 완료"; return 0; }
  mkdir -p $O
  say "$TAG 시작"
  $PY ./l1_train_v10.py "${common[@]}" --out $O "$@" > $O/run.log 2>&1
  local RC=$?
  say "$TAG rc=$RC $( [ -s $O/L1_DONE ] && echo L1_DONE || tail -2 $O/run.log | tr '\n' ' ' | cut -c1-140 )"
}

case "$CH" in
  *seventeen*) run sens_phi17_f0 --phi-columns v9-seventeen --c1-ecdf-reference training ;;
  *) say "17열 선택지 없음 — S1 건너뜀 (선택지: $CH)" ;;
esac
run sens_ecdfeval_f0 --phi-columns v9-eighteen --c1-ecdf-reference evaluation

say "SENS_NIGHT_DONE"
echo ok > $R/sens_night.done
