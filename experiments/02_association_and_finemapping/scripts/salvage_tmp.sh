#!/usr/bin/env bash
: "${CONDA_PREFIX:?Set CONDA_PREFIX}"
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
export PATH=${CONDA_PREFIX}/bin:$PATH
W=${PROJECT_ROOT}/work/run_ourfm
for T in $W/regcache/*.vcf.gz.tmp; do
  V=${T%.tmp}
  EOF_=$(tail -c 28 $T | od -An -tx1 | tr -d ' \n')
  [ "$EOF_" = "1f8b08040000000000ff0600424302001b0003000000000000000000" ] || { echo "SKIP no-EOF $(basename $T)"; continue; }
  mv $T $V; rm -f $T.csi
  bcftools index -c -f $V || { echo "INDEX FAIL $(basename $V)"; continue; }
  N=$(bcftools index -n $V); [ "${N:-0}" -gt 0 ] && echo "ok $N" > $V.done
  echo "salvaged n=$N"
done
