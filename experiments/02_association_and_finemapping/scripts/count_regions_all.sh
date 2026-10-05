#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
W=${PROJECT_ROOT}/work/run_ourfm
for N in $(seq 1 22); do for T in tchl htn dm lip; do
  if [ -s $W/regions/$T.chr$N.regions.tsv ]; then n=$(( $(wc -l < $W/regions/$T.chr$N.regions.tsv) - 1 )); echo "[regions $T chr$N] existing n_region=$n"; else python3 $W/scripts/fm_regions.py $T $N; fi
done; done
echo COUNT_DONE
