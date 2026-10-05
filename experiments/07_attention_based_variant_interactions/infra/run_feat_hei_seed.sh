#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT to the working data root}"
[[ "$PROJECT_ROOT" =~ [^[:space:]] ]] || { printf "%s\n" "PROJECT_ROOT must be nonblank" >&2; exit 2; }
CODE_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
export CODE_ROOT
set -uo pipefail
cd ${PROJECT_ROOT}/work/prs; D=out/hei_seed/feat; mkdir -p $D
[ -s out/hei/sel_seed/select_seed.done ] || { echo "no select_seed.done"; exit 1; }
if [ ! -s $D/token_variants.tsv ]; then
  awk -F'\t' 'NR>1 && !s[$3]++ {print $3"\t"$2"\t"$4"\t"$1}' out/hei/sel_seed/gene_tokens.tsv | (echo -e "key\tchr\tpos38\tdomain"; cat) > $D/token_variants.tsv.tmp && mv $D/token_variants.tsv.tmp $D/token_variants.tsv
fi
echo "variants $(($(wc -l < $D/token_variants.tsv)-1))"
cut -f2 $D/token_variants.tsv | tail -n+2 | sort -u | xargs -P 6 -I{} bash -c '[ -s out/hei_seed/feat/feat_chr{}.tsv.gz ] || "${PYTHON_BIN:-python3}" "${CODE_ROOT}/data_prep/f_features_hei_seed.py" {} out/hei_seed/feat/token_variants.tsv >> logs/feat_hei_seed.jsonl 2>> logs/feat_hei_seed.err'
n=$(ls $D/feat_chr*.tsv.gz 2>/dev/null | wc -l); echo "feat $n $(date -Is)"; [ $n -eq 22 ] && echo "ok $(date -Is)" > $D/feat.done
