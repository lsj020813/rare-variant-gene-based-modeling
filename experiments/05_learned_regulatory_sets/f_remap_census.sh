#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
W=${PROJECT_ROOT}/work/fset; mkdir -p $W/out $W/logs
zcat ${PROJECT_ROOT}/work/ref/remap/remap2022_nr_macs2_hg38_v1_0.bed.gz \
 | awk -F'\t' '{split($4,a,":"); n[a[1]]++} END{for(t in n) print n[t]"\t"t}' \
 | sort -k1,1nr > $W/out/remap_tf_counts.tsv.tmp
mv $W/out/remap_tf_counts.tsv.tmp $W/out/remap_tf_counts.tsv
head -60 $W/out/remap_tf_counts.tsv | cut -f2 > $W/out/remap_top_tfs.txt
wc -l < $W/out/remap_tf_counts.tsv > $W/out/remap_tf_total.txt
echo "remap census done $(date -Is) tfs=$(cat $W/out/remap_tf_total.txt)" > $W/out/remap_census.done
