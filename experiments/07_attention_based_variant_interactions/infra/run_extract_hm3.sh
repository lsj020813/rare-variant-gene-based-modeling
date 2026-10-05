#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT to the working data root}"
[[ "$PROJECT_ROOT" =~ [^[:space:]] ]] || { printf "%s\n" "PROJECT_ROOT must be nonblank" >&2; exit 2; }
CODE_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
export CODE_ROOT
cd ${PROJECT_ROOT}/work/prs
seq 1 22 | xargs -P 8 -I{} bash "${CODE_ROOT}/baselines/prs_cs/extract_hm3.sh" {}
echo "ALL $(date -Is)" > geno_hm3/all.done
