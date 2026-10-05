#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT to the working data root}"
[[ "$PROJECT_ROOT" =~ [^[:space:]] ]] || { printf "%s\n" "PROJECT_ROOT must be nonblank" >&2; exit 2; }
CODE_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
export CODE_ROOT
cd ${PROJECT_ROOT}/work/prs
cut -f4 out/loci.bed | xargs -P 8 -I{} bash -c '[ -s private/loci/{}.npz ] && exit 0; "${PYTHON_BIN:-python3}" "${CODE_ROOT}/data_prep/build_locus.py" {} >> logs/build_loci.jsonl 2>> logs/build_loci.err'
echo "ALL $(date -Is) $(ls private/loci/*.npz | wc -l)" > out/build_loci.done
