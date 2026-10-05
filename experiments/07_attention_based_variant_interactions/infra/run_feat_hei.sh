#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT to the working data root}"
[[ "$PROJECT_ROOT" =~ [^[:space:]] ]] || { printf "%s\n" "PROJECT_ROOT must be nonblank" >&2; exit 2; }
CODE_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
export CODE_ROOT
cd ${PROJECT_ROOT}/work/prs
cut -f2 out/hei/feat/token_variants.tsv | tail -n+2 | sort -u | xargs -P 6 -I{} bash -c '[ -s out/hei/feat/feat_chr{}.tsv.gz ] || "${PYTHON_BIN:-python3}" "${CODE_ROOT}/data_prep/f_features_hei.py" {} out/hei/feat/token_variants.tsv >> logs/feat_hei.jsonl 2>> logs/feat_hei.err'
n=$(ls out/hei/feat/feat_chr*.tsv.gz 2>/dev/null | wc -l); echo "feat $n $(date -Is)"; [ $n -eq 22 ] && echo ok > out/hei/feat/feat.done
