#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT to the working data root}"
[[ "$PROJECT_ROOT" =~ [^[:space:]] ]] || { printf "%s\n" "PROJECT_ROOT must be nonblank" >&2; exit 2; }
CODE_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
export CODE_ROOT
cd ${PROJECT_ROOT}/work/prs; PY="${PYTHON_BIN:-python3}"
while [ ! -s out/hei/sel/select.done ]; do sleep 30; done
[ -s out/hei/sel/selection_summary.json ] || { echo "selection failed"; exit 1; }
ls out/hei/sel/extract_chr*.keys | sed 's/.*extract_chr//; s/.keys//' | sort -n | xargs -P 6 -I{} $PY "${CODE_ROOT}/data_prep/hei_extract.py" {}
n=$(ls out/hei/meta/chr*.done | wc -l); m=$(ls out/hei/sel/extract_chr*.keys | wc -l); echo "extract $n/$m $(date -Is)"; [ $n -eq $m ] && echo ok > out/hei/extract.done
