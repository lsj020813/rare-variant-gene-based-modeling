#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT to the working data root}"
[[ "$PROJECT_ROOT" =~ [^[:space:]] ]] || { printf "%s\n" "PROJECT_ROOT must be nonblank" >&2; exit 2; }
CODE_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
export CODE_ROOT
set -uo pipefail
cd ${PROJECT_ROOT}/work/prs; PY="${PYTHON_BIN:-python3}"; O=out/hei_seed
[ -s out/hei/sel_seed/select_seed.done ] || { echo "no select_seed.done"; exit 1; }
mkdir -p $O
if [ ! -s $O/extract.done ]; then
  for kf in out/hei/sel_seed/extract_chr*.keys; do N=$(basename $kf .keys | sed 's/extract_chr//'); nk=$(grep -c . $kf); for k in $(seq 0 $(( (nk+4999)/5000 - 1 ))); do echo "$N $k $nk"; done; done | sort -k3,3nr | awk '{print $1" "$2}' > $O/extract_jobs.txt
  echo "dosage jobs $(wc -l < $O/extract_jobs.txt) $(date -Is)"
  cat $O/extract_jobs.txt | xargs -P 20 -n 2 $PY "${CODE_ROOT}/data_prep/hei_seed_extract_part.py" >> logs/hei_seed_parts.log 2>> logs/hei_seed_parts.err
  $PY "${CODE_ROOT}/data_prep/hei_seed_extract_merge.py" > logs/hei_seed_merge.log 2>&1 || { echo "dosage merge failed"; exit 2; }
  n=$(ls $O/meta/chr*.done | wc -l); m=$(ls out/hei/sel_seed/extract_chr*.keys | wc -l); echo "dosage $n/$m $(date -Is)"
  [ $n -eq $m ] && echo "ok $(date -Is)" > $O/extract.done || exit 3
fi
mkdir -p $O/gt
if [ ! -s $O/gt/extract.done ]; then
  ls $O/gt/keys_chr*.txt >/dev/null 2>&1 || awk -F'\t' -v d=$O/gt 'NR>1 && $6<128 {print $3 > (d"/keys_chr"$2".txt.tmp")}' out/hei/sel_seed/gene_tokens.tsv
  for f in $O/gt/keys_chr*.txt.tmp; do [ -e "$f" ] || continue; sort -u $f | sort -t: -k2,2n -k3,3 -k4,4 > ${f%.tmp}; rm $f; done
  echo "gt keys $(cat $O/gt/keys_chr*.txt | wc -l)"
  for kf in $O/gt/keys_chr*.txt; do N=$(basename $kf .txt | sed 's/keys_chr//'); nk=$(grep -c . $kf); for k in $(seq 0 $(( (nk+4999)/5000 - 1 ))); do echo "$N $k $nk"; done; done | sort -k3,3nr | awk '{print $1" "$2}' > $O/gt/jobs.txt
  echo "gt jobs $(wc -l < $O/gt/jobs.txt) $(date -Is)"
  cat $O/gt/jobs.txt | xargs -P 20 -n 2 $PY "${CODE_ROOT}/data_prep/hei_seed_extract_gt_part.py" >> logs/hei_seed_gt_parts.log 2>> logs/hei_seed_gt_parts.err
  $PY "${CODE_ROOT}/data_prep/hei_seed_extract_gt_merge.py" > logs/hei_seed_gt_merge.log 2>&1 || { echo "gt merge failed"; exit 4; }
  n=$(ls $O/gt/chr*.done | wc -l); m=$(ls $O/gt/keys_chr*.txt | wc -l); echo "gt $n/$m $(date -Is)"
  [ $n -eq $m ] && echo "ok $(date -Is)" > $O/gt/extract.done || exit 5
fi
echo "ALL DONE $(date -Is)"
