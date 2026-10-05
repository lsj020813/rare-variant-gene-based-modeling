: "${PROJECT_ROOT:?Set PROJECT_ROOT to the working data root}"
[[ "$PROJECT_ROOT" =~ [^[:space:]] ]] || { printf "%s\n" "PROJECT_ROOT must be nonblank" >&2; exit 2; }
CODE_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
export CODE_ROOT
cd ${PROJECT_ROOT}/work/prs; PY="${PYTHON_BIN:-python3}"
for T in $THREADS; do OMP_NUM_THREADS=$T MKL_NUM_THREADS=$T timeout 900 nice -n 10 $PY "${CODE_ROOT}/models/attention/bench_cpu_gpu2.py" cpu $T 2>&1 | tail -1; echo "load $(cut -d' ' -f1 /proc/loadavg) $(date +%H:%M:%S)"; done
echo SCALE_END
