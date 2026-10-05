#!/bin/bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -u
ATT=${PROJECT_ROOT}/work/fset/att
OUT=$ATT/hqcache
PY=python3
mkdir -p $OUT
cd $ATT
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export TMPDIR=${PROJECT_ROOT}/work/tmp
export R2_MIN=${R2_MIN:-0.9}
N=${NWORK:-10}
echo "CACHE_START R2_MIN=$R2_MIN N=$N $(date -Is)" | tee -a $OUT/run.log
for ch in $(seq 1 22); do
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 5; done
  ( nice -n 19 ionice -c3 $PY cache_hq.py $ch > $OUT/chr${ch}.log 2>&1 ; \
    echo "chr${ch} exit=$?" >> $OUT/exit.log ) &
done
wait
echo "CACHE_DONE $(date -Is)" | tee -a $OUT/run.log
$PY - <<'PYEOF' >> $OUT/run.log 2>&1
import os as _env_os, re as _env_re
def _env_path(value):
    return _env_re.sub(r"\$\{([A-Z][A-Z0-9_]*)\}", lambda match: _env_os.environ[match.group(1)], value)
import glob, json
tot = 0; gb = 0.0; bad = []
for f in sorted(glob.glob(_env_path('${PROJECT_ROOT}/work/fset/att/hqcache/hq_chr*.json'))):
    j = json.load(open(f)); tot += j['n_matched']; gb += j['gb']
    if j['match_rate'] < 0.95: bad.append((j['chrom'], j['match_rate']))
print('총 캐시 변이 %d개, %.1f GB, 매칭률 95%% 미만 염색체: %s' % (tot, gb, bad))
PYEOF
touch $OUT/cache.done
