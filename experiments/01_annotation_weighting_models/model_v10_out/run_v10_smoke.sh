#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
export TMPDIR="$PWD/.scratch" TMP="$PWD/.scratch" TEMP="$PWD/.scratch"
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
mkdir -p -- "$TMPDIR"
ulimit -c 0
ulimit -v 8388608
v10_python=python3
v10_fixture=${1:-.smoke/fixture_v10}
v10_output=${2:-.smoke/verified}
if [[ ! -f "$v10_fixture/FIXTURE.json" ]]; then
  "$v10_python" -B l1_train_v10.py --make-smoke-fixture "$v10_fixture" \
    --out "$v10_output/setup" --threads 4 --memory-gb 8 --tmpdir "$TMPDIR"
fi
exec "$v10_python" -B l1_train_v10.py --smoke-suite "$v10_fixture" \
  --out "$v10_output" --threads 4 --memory-gb 8 --tmpdir "$TMPDIR"
