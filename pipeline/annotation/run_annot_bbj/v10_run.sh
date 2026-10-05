#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -u
NAME="${1:?name}"; SUF="${2-}"; shift 2
cd ${PROJECT_ROOT}/work/run_annot_bbj
OUT=${PROJECT_ROOT}/work/run_band15/model_v10_out/e6_runs/$NAME; mkdir -p "$OUT"
L1=${PROJECT_ROOT}/work/ref/annot/bbj_l1
PY=python3
CODE=${PROJECT_ROOT}/work/run_band15/model_v10_out/l1_train_v10.py
export TMPDIR=${PROJECT_ROOT}/work/run_band15/model_v10_out/.scratch_e6/$NAME; export TMP=$TMPDIR TEMP=$TMPDIR; mkdir -p "$TMPDIR"
export PYTHONDONTWRITEBYTECODE=1
echo "$(date '+%F %T') START $NAME args: $*" >> "$OUT/run.meta"
( "$PY" -B "$CODE" --out "$OUT" --bbj-annot "$L1/bbj_annot$SUF.tsv.gz" --bbj-pip "$L1/bbj_pip$SUF" --traits "$L1/traits$SUF.tsv" \
    --annot-stats ${PROJECT_ROOT}/work/ref/annot/cache/fm_all.stats.json \
    --threads 4 --memory-gb 60 --tmpdir "$TMPDIR" --allow-concurrent-heavy \
    --c1-columns v9-eight --unseen-intercept-rule mean-trained --export-phi --export-intercept-rule mean-trained \
    --alpha 0.05 --power-target 0.8 --type1-max 0.10 "$@" > "$OUT/run.log" 2>&1; echo "rc=$?" > "$OUT/rc.txt"; echo "$(date '+%F %T') END rc=$(cat $OUT/rc.txt)" >> "$OUT/run.meta"
  "$PY" - "$OUT/RUN.json" <<'PYEOF'
import json, sys, os
p = sys.argv[1]
d = json.load(open(p)) if os.path.exists(p) else {}
d["e6_policy"] = {"ram_cap_gb_per_job": 60, "ram_cap_basis": "relayed by parent monitor frame 2026-09-09 07:4x as a user decision (38 = core cap, not RAM); direct user confirmation pending (STAR)", "core_cap": 38,
                  "watchdog": "wd_rss.sh actual-RSS, 2 consecutive samples > cap -> kill; no trend extrapolation", "rlimit_as_gb": 60}
json.dump(d, open(p, "w"), indent=1)
PYEOF
) &
sleep 5
bash ${PROJECT_ROOT}/work/run_annot_bbj/wd_rss.sh "l1_train_v1[0].py.*e6_runs/$NAME " "$OUT/wd_rss.log" 60 6 &
wait
