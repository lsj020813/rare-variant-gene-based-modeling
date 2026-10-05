#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
cd ${PROJECT_ROOT}/work/phi_gate; PY=python3
while pgrep -f 'python final_phi_gate_v9[.]py chr ' >/dev/null; do sleep 60; done
for N in 12 2; do rc=$(grep -o 'exit=[0-9]*' logs/v9_chr$N.log | tail -1 | cut -d= -f2); echo "chr $N rc=$rc $(date -Is)" | tee -a logs/chain_v9.log
  [ "$rc" = "0" ] || { echo "CHAIN STOPPED: chr $N failed" | tee -a logs/chain_v9.log; exit 1; }; done
echo "START pool $(date -Is)" | tee -a logs/chain_v9.log
$PY final_phi_gate_v9.py pool 2 12 > logs/v9_pool.log 2>&1; rc=$?
echo "END pool rc=$rc $(date -Is)" | tee -a logs/chain_v9.log
[ $rc -eq 0 ] && echo "chain ok $(date -Is)" > out/chain_v9.done
