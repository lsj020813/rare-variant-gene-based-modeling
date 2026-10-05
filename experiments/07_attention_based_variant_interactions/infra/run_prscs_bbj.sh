#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT to the working data root}"
[[ "$PROJECT_ROOT" =~ [^[:space:]] ]] || { printf "%s\n" "PROJECT_ROOT must be nonblank" >&2; exit 2; }
CODE_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
export CODE_ROOT
cd ${PROJECT_ROOT}/work/prs
for N in $(seq 1 22); do while [ ! -s geno_hm3/chr$N.done ]; do sleep 30; done; done
seq 1 22 | xargs -P 11 -I{} bash "${CODE_ROOT}/baselines/prs_cs/prscs_chr.sh" {} bbj ref/ldblk_1kg_eas out/bbj_tc_prscs.txt 156719
ls prscs_bbj/*.done | wc -l > prscs_bbj/all.count; echo "ALL $(date -Is)" >> logs/prscs_bbj.log
