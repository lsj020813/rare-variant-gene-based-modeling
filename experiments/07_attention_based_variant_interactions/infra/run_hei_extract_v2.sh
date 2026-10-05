#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT to the working data root}"
[[ "$PROJECT_ROOT" =~ [^[:space:]] ]] || { printf "%s\n" "PROJECT_ROOT must be nonblank" >&2; exit 2; }
CODE_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
export CODE_ROOT
cd ${PROJECT_ROOT}/work/prs; PY="${PYTHON_BIN:-python3}"
for kf in out/hei/sel/extract_chr*.keys; do N=$(basename $kf .keys | sed 's/extract_chr//'); nk=$(grep -c . $kf); for k in $(seq 0 $(( (nk+4999)/5000 - 1 ))); do echo "$N $k $nk"; done; done | sort -k3,3nr | awk '{print $1" "$2}' > out/hei/extract_jobs.txt
echo "jobs $(wc -l < out/hei/extract_jobs.txt)"
cat out/hei/extract_jobs.txt | xargs -P 20 -n 2 $PY "${CODE_ROOT}/data_prep/hei_extract_part.py" >> logs/hei_extract_parts.log 2>> logs/hei_extract_parts.err
$PY "${CODE_ROOT}/data_prep/hei_extract_merge.py" > logs/hei_extract_merge.log 2>&1 || exit 2
n=$(ls out/hei/meta/chr*.done | wc -l); m=$(ls out/hei/sel/extract_chr*.keys | wc -l); echo "extract $n/$m $(date -Is)"; [ $n -eq $m ] && echo ok > out/hei/extract.done
