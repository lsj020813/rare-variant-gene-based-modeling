#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
W=${PROJECT_ROOT}/work/gate1; PY=python3
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "CHAIN START $(date -Is) pid=$$"
for i in $(seq 1 240); do [ -s $W/out/windows.done ] && break; sleep 30; done
[ -s $W/out/windows.done ] || { echo "ABORT windows.done missing"; exit 1; }
echo "windows ready $(date -Is): $(cat $W/out/windows.done)"
cds_one(){ N=$1; W=${PROJECT_ROOT}/work/gate1; PY=python3
  [ -s $W/out/vpos_chr$N.tsv.gz ] && { echo "HAVE vpos chr$N"; return 0; }
  nice -n 19 ionice -c3 $PY $W/gate1_cdsflag.py $N > $W/logs/cds_chr$N.json 2> $W/logs/cds_chr$N.err
  echo "cds chr$N rc=$? $(head -c 200 $W/logs/cds_chr$N.json)"; }
export -f cds_one
seq 22 -1 1 | xargs -P 6 -I{} bash -c 'cds_one {}'
echo "cdsflag done $(date -Is) files=$(ls $W/out/vpos_chr*.tsv.gz 2>/dev/null | wc -l)"
$PY $W/gate1_sample_genes.py > $W/logs/gene_sample.json 2> $W/logs/gene_sample.err
echo "sample: $(cat $W/logs/gene_sample.json)"
[ -s $W/out/gene_sample.csv ] || { echo "ABORT gene_sample.csv missing"; head -5 $W/logs/gene_sample.err; exit 1; }
ctx_one(){ N=$1; W=${PROJECT_ROOT}/work/gate1; PY=python3
  [ -s $W/out/genes_chr$N.txt ] || return 0
  [ -s $W/out/ctx_chr$N.json ] && { echo "HAVE ctx chr$N"; return 0; }
  nice -n 19 ionice -c3 $PY $W/gate1_context_v2.py $N $W/out/genes_chr$N.txt > $W/logs/ctx_chr$N.log 2> $W/logs/ctx_chr$N.err
  echo "ctx chr$N rc=$? $(tail -1 $W/logs/ctx_chr$N.log)"; }
export -f ctx_one
seq 22 -1 1 | xargs -P 12 -I{} bash -c 'ctx_one {}'
n=$(ls $W/out/ctx_chr*.json 2>/dev/null | wc -l)
echo "CHAIN END $(date -Is) ctx_files=$n"
echo "gate1 chain complete $(date -Is) ctx=$n" > $W/out/gate1_chain.done
