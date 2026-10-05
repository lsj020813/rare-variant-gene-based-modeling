#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
cd ${PROJECT_ROOT}/work/run_band15/model_v10_out
export TMPDIR="$PWD/.scratch" TMP="$PWD/.scratch" TEMP="$PWD/.scratch"
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
mkdir -p "$TMPDIR"; ulimit -c 0; ulimit -v 8388608
PY=python3
F=.smoke/fixture_e6; O=.smoke/verified_e6
rm -rf "$F" "$O"
"$PY" -B l1_train_v10.py --make-smoke-fixture "$F" --out "$O/setup" --threads 4 --memory-gb 8 --tmpdir "$TMPDIR" --allow-concurrent-heavy || exit 10
"$PY" -B l1_train_v10.py --smoke-suite "$F" --out "$O" --threads 4 --memory-gb 8 --tmpdir "$TMPDIR" --allow-concurrent-heavy
