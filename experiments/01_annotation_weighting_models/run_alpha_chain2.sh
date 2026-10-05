#!/usr/bin/env bash
: "${CONDA_INIT_SCRIPT:?Set CONDA_INIT_SCRIPT}"
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -o pipefail
B=${PROJECT_ROOT}/work/run_l3b; R=${PROJECT_ROOT}/work/ref
source ${CONDA_INIT_SCRIPT} && conda activate "${CONDA_ENV:?Set CONDA_ENV}"
export TMPDIR=${PROJECT_ROOT}/work/tmp TMP=${PROJECT_ROOT}/work/tmp TEMP=${PROJECT_ROOT}/work/tmp
PY="nice -n 19 ionice -c3 python3"
cd $B
for A in 0.1 0.25 0.5; do
  TAG=alpha${A/./p}
  echo "=== build $TAG $(date '+%m-%d %H:%M:%S') ==="
  if [ ! -s out/gf_${TAG}_manifest.done ]; then
    $PY mk_gf_alpha.py --alpha $A --g-phi out/G_phi.txt --maf-tsv out/maf19.tsv --out out/G_$TAG.txt --threads 4 || { echo "CHAIN ABORTED build $TAG"; exit 2; }
  fi
  ng=$(awk '$2=="var"{print $1}' out/G_$TAG.txt | sort -u | wc -l); ng0=$(awk '$2=="var"{print $1}' out/G_phi.txt | sort -u | wc -l)
  ns=$(awk '$2=="weight"{n+=NF-2} END{print n}' out/G_$TAG.txt);   ns0=$(awk '$2=="weight"{n+=NF-2} END{print n}' out/G_phi.txt)
  [ "$ng" = "$ng0" ] && [ "$ns" = "$ns0" ] || { echo "GATE FAIL $TAG genes $ng/$ng0 slots $ns/$ns0"; exit 5; }
  echo "[$TAG] gate ok genes=$ng slots=$ns"
  if [ "$(ls out/chunks/G_$TAG.part*.txt 2>/dev/null | wc -l)" -eq 0 ]; then
    $PY split_chunks.py --group-files out/G_$TAG.txt --genes-per-chunk 3 --out out/chunks || { echo "CHAIN ABORTED split $TAG"; exit 2; }
  fi
  nc=$(ls out/chunks/G_$TAG.part*.txt | wc -l); nc0=$(ls out/chunks/G_phibeta.part*.txt | wc -l)
  [ "$nc" = "$nc0" ] || { echo "GATE FAIL $TAG chunks $nc/$nc0"; exit 5; }
done
for A in 0.1 0.25 0.5; do
  TAG=alpha${A/./p}
  [ -s out/$TAG.merged.done ] && { echo "=== $TAG already merged, skip"; continue; }
  echo "=== $TAG dispatch $(date '+%m-%d %H:%M:%S') load=$(cut -d' ' -f1 /proc/loadavg) ==="
  rm -f out/STOP
  bash run_arm_seq2.sh $TAG 16 45 70; rc=$?
  echo "=== $TAG rc=$rc $(date '+%m-%d %H:%M:%S') ==="
  [ $rc -ne 0 ] && { echo "CHAIN ABORTED at $TAG"; exit $rc; }
done
echo "ALPHA_CHAIN_DONE $(date '+%m-%d %H:%M:%S')"
