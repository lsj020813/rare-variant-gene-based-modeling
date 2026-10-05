#!/bin/bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -o pipefail
export TMPDIR=${PROJECT_ROOT}/work/tmp TMP=${PROJECT_ROOT}/work/tmp TEMP=${PROJECT_ROOT}/work/tmp
W=${PROJECT_ROOT}/work/run_meth; mkdir -p $W/out $W/logs
V=${HEXA_VCF_FILE:?Set HEXA_VCF_FILE}
echo "START $(date -Is) host=$(hostname) vcf_size=$(stat -c %s $V)"
nice -n 19 ionice -c3 python $W/gate1_carriers.py $V $W/out/gate1_chr22.json
rc=$?
echo "END $(date -Is) rc=$rc"
[ $rc -eq 0 ] && echo "gate1 chr22 rc=0 $(date -Is)" > $W/out/gate1_chr22.done
