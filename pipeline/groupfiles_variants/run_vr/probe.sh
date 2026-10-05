#!/usr/bin/env bash
set -uo pipefail
BCF=${BCFTOOLS:-bcftools}
SRC=${GENOTYPE_DIR:?Set GENOTYPE_DIR}/${GENOTYPE_FILENAME_PREFIX:?Set GENOTYPE_FILENAME_PREFIX}22${GENOTYPE_FILENAME_SUFFIX:?Set GENOTYPE_FILENAME_SUFFIX}
ls -l $SRC | awk '{printf "source %.1fGB\n", $5/1e9}'
echo "--- INFO fields present ---"
$BCF view -h $SRC 2>/dev/null | grep -E '^##INFO' | head -6
echo "--- first 200k records: MAF distribution ---"
$BCF query -f '%INFO/MAF\t%INFO/R2\n' $SRC 2>/dev/null | head -200000 | awk -F'\t' '{m=$1+0; r=$2+0; n++;
  if(m>5.72e-05 && m<=1.1724e-04){a++; if(r>=0.3) ar++}
  else if(m<=5.72e-05) z++;
  else b++}
 END{printf "  scanned %d\n  MAC<=10 : %d\n  MAC 10-20.5: %d  (R2>=0.3: %d)\n  MAC>20.5: %d\n", n,z,a,ar,b}'
echo "SMOKE_DONE"
