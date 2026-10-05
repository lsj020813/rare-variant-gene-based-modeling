#!/bin/bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -o pipefail
export TMPDIR=${PROJECT_ROOT}/work/tmp
cd ${PROJECT_ROOT}/work/run_trackB/ref
curl -sL -o chr19.hg19.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr19.fa.gz && md5sum chr19.hg19.fa.gz && ls -la chr19.hg19.fa.gz
curl -sL -o hg19.md5sum.txt https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/md5sum.txt; grep "chr19.fa.gz" hg19.md5sum.txt
gunzip -kf chr19.hg19.fa.gz && ls -la chr19.hg19.fa
curl -sL -o targets_human.txt https://raw.githubusercontent.com/calico/basenji/master/manuscripts/cross2020/targets_human.txt && wc -l targets_human.txt && md5sum targets_human.txt && head -2 targets_human.txt
echo EXIT $?; date; touch fetch_ref.done
