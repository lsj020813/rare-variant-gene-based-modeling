#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT to the working data root}"
[[ "$PROJECT_ROOT" =~ [^[:space:]] ]] || { printf "%s\n" "PROJECT_ROOT must be nonblank" >&2; exit 2; }
CODE_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
export CODE_ROOT
cd ${PROJECT_ROOT}/work/prs; PY="${PYTHON_BIN:-python3}"
mkdir -p out/hei/gt
awk -F'\t' 'NR>1 && $6<128 {print $2"\t"$3}' out/hei/sel/gene_tokens.tsv | sort -u | awk -F'\t' '{print $2 > "out/hei/gt/keys_chr"$1".txt"}'
for f in out/hei/gt/keys_chr*.txt; do sort -u $f | sort -t: -k2,2n -o $f; done
echo "keys $(cat out/hei/gt/keys_chr*.txt | wc -l)"
for kf in out/hei/gt/keys_chr*.txt; do N=$(basename $kf .txt | sed 's/keys_chr//'); nk=$(grep -c . $kf); for k in $(seq 0 $(( (nk+4999)/5000 - 1 ))); do echo "$N $k $nk"; done; done | sort -k3,3nr | awk '{print $1" "$2}' > out/hei/gt/jobs.txt
echo "jobs $(wc -l < out/hei/gt/jobs.txt)"
cat out/hei/gt/jobs.txt | xargs -P 20 -n 2 $PY "${CODE_ROOT}/data_prep/hei_extract_gt_part.py" >> logs/hei_gt_parts.log 2>> logs/hei_gt_parts.err
$PY "${CODE_ROOT}/data_prep/hei_extract_gt_merge.py" > logs/hei_gt_merge.log 2>&1 || exit 2
n=$(ls out/hei/gt/chr*.done | wc -l); m=$(ls out/hei/gt/keys_chr*.txt | wc -l); echo "gt $n/$m $(date -Is)"; [ $n -eq $m ] && echo ok > out/hei/gt/extract.done
