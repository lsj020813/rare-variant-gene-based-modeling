#!/usr/bin/env bash
for N in $(seq 1 22); do python3 assemble.py $N > fm_$N.log 2>&1 || echo "[chr$N] FAIL"; done
echo FM2_DONE
