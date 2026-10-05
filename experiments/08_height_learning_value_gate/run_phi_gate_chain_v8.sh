#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
cd ${PROJECT_ROOT}/work/phi_gate; PY=python3
while pgrep -f 'python final_phi_gate_v8[.]py chr 12' >/dev/null; do sleep 60; done
rc12=$(grep -o 'exit=[0-9]*' logs/chr12_test.log | tail -1 | cut -d= -f2)
echo "chr 12 finished rc=$rc12 $(date -Is)" | tee -a logs/chain_v8.log
[ "$rc12" = "0" ] || { echo "CHAIN STOPPED: chr 12 failed" | tee -a logs/chain_v8.log; exit 1; }
for step in "chr 2" "pool 2 12"; do
  tag=$(echo $step | tr ' ' '_'); L=logs/gate_v8_${tag}.log
  echo "START $step $(date -Is)" | tee -a logs/chain_v8.log
  $PY final_phi_gate_v8.py $step > $L 2>&1; rc=$?
  echo "END $step rc=$rc $(date -Is)" | tee -a logs/chain_v8.log
  [ $rc -eq 0 ] || { echo "CHAIN STOPPED at $step" | tee -a logs/chain_v8.log; exit $rc; }
done
echo "chain ok $(date -Is)" > out/chain_v8.done
