#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT to the working data root}"
[[ "$PROJECT_ROOT" =~ [^[:space:]] ]] || { printf "%s\n" "PROJECT_ROOT must be nonblank" >&2; exit 2; }
CODE_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
export CODE_ROOT
cd ${PROJECT_ROOT}/work/prs
( echo -e "key\tchr\tpos38\tdomain"; for f in out/loci_meta/L*.tsv; do L=$(basename $f .tsv); awk -v L=$L 'NR>1{print $1"\t"$2"\t"$4"\t"L}' $f; done ) > out/feat/locus_variants.tsv
cut -f2 out/feat/locus_variants.tsv | tail -n+2 | sort -u | xargs -P 8 -I{} bash -c '[ -s out/feat/feat_chr{}.tsv.gz ] || "${PYTHON_BIN:-python3}" "${CODE_ROOT}/data_prep/f_features_prs.py" {} out/feat/locus_variants.tsv >> logs/feat.jsonl 2>> logs/feat.err'
echo "ALL $(date -Is)" > out/feat/feat.done
