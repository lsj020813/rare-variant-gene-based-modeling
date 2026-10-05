#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
W=${PROJECT_ROOT}/work/fset; D=$W/out/remap_by_chr; mkdir -p $D
[ -s $D/chr22.bed ] && { echo "HAVE remap split"; exit 0; }
zcat ${PROJECT_ROOT}/work/ref/remap/remap2022_nr_macs2_hg38_v1_0.bed.gz \
 | awk -F'\t' -v d="$D" '$1 ~ /^chr([0-9]+)$/ {print > (d"/"$1".bed.tmp")}'
for f in $D/*.bed.tmp; do mv "$f" "${f%.tmp}"; done
echo "remap split done $(date -Is) files=$(ls $D | wc -l)" > $W/out/remap_split.done
