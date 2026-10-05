#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
cd ${PROJECT_ROOT}/work/phi_gate; PY=python3
for step in "chr 12" "chr 2" "pool 2 12"; do
  tag=$(echo $step | tr ' ' '_'); L=logs/gate_${tag}.log
  echo "START $step $(date -Is)" | tee -a logs/chain.log
  $PY final_phi_gate_v7.py $step > $L 2>&1; rc=$?
  echo "END $step rc=$rc $(date -Is)" | tee -a logs/chain.log
  [ $rc -eq 0 ] || { echo "CHAIN STOPPED at $step" | tee -a logs/chain.log; exit $rc; }
done
echo "chain ok $(date -Is)" > out/chain.done
