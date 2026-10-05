#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
cd ${PROJECT_ROOT}/work/run_band15
while true; do n=$(ls ${PROJECT_ROOT}/work/ref15/annot/ds/*.done 2>/dev/null | wc -l); pool=$(pgrep -fc 'ds_stor[e]' || true)
  if [ "$n" -eq 22 ] && [ "${pool:-0}" -eq 0 ]; then break; fi; sleep 60; done
echo "[waiter] $(date '+%m-%d %H:%M') DS 22/22, pool idle — checking totals"
python3 - <<'EOF'
import os as _env_os, re as _env_re
def _env_path(value):
    return _env_re.sub(r"\$\{([A-Z][A-Z0-9_]*)\}", lambda match: _env_os.environ[match.group(1)], value)
import json, glob
t = 0
for f in glob.glob(_env_path("${PROJECT_ROOT}/work/ref15/annot/ds/chr*.ds.npz.stats.json")):
    d = json.load(open(f)); t += d["bytes"]
print(f"[waiter] DS total {t/1e9:.1f} GB across {len(glob.glob(_env_path('${PROJECT_ROOT}/work/ref15/annot/ds/chr*.ds.npz.stats.json')))} files")
EOF
echo "[waiter] launching l1_chain_v7 (band15)"
setsid nohup bash l1_chain_v7.sh > l1_chain_v7.log 2>&1 </dev/null &
echo WAITER_LAUNCHED
