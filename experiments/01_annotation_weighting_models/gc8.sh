#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
cd ${PROJECT_ROOT}/work/run_band15
for arm in spline linear nn cluster; do
  bash rt.sh gc8_$arm.log --smoke 22 --arm $arm --gradcheck-only --threads 6 --known-genes ${PROJECT_ROOT}/work/ref15/annot/l1/known_genes_truthB.tsv
  echo "=== $arm exit $?"
done
echo GC8_ALL_DONE
