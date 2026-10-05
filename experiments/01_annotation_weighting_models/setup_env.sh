#!/bin/bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -o pipefail
export TMPDIR=${PROJECT_ROOT}/work/tmp TMP=${PROJECT_ROOT}/work/tmp TEMP=${PROJECT_ROOT}/work/tmp
cd ${PROJECT_ROOT}/work/run_trackB
PY=python3
$PY -m venv --system-site-packages venv || exit 1
./venv/bin/pip install --no-cache-dir "einops>=0.6" pyfaidx polars "transformers<4.47" huggingface_hub 2>&1 | tail -5
./venv/bin/pip install --no-cache-dir --no-deps enformer-pytorch 2>&1 | tail -3
./venv/bin/python -c "import enformer_pytorch, transformers, einops, torch; print('OK', enformer_pytorch.__version__ if hasattr(enformer_pytorch,'__version__') else 'v?', transformers.__version__, einops.__version__, torch.__version__, torch.cuda.is_available())"
echo EXIT $? ; date
touch setup_env.done
