#!/usr/bin/env bash
set -euo pipefail
: "${PROJECT_ROOT:?Set PROJECT_ROOT to a generic workspace root}"
OUT="${PROJECT_ROOT}/work/ref/saige_step1"
GRM=${GRM_FILE:?Set GRM_FILE}
GRM_IDS=${GRM_SAMPLE_IDS_FILE:?Set GRM_SAMPLE_IDS_FILE}
PLINK=${PLINK_PREFIX:?Set PLINK_PREFIX}
mkdir -p "$OUT"
SAIGE_RUNNER=${SAIGE_RUNNER:-Rscript}
SAIGE_STEP1_SCRIPT=${SAIGE_STEP1_SCRIPT:?Set SAIGE_STEP1_SCRIPT}
for T in htn dm lip; do
  P="${PHENO_V2_DIR:?Set PHENO_V2_DIR}/${T}_v2.tsv"
  test -s "$P"
  if [ -s "$OUT/${T}_v2.rda" ] && [ -s "$OUT/${T}_v2.varianceRatio.txt" ]; then echo "[${T}] done, skip"; continue; fi
  echo "[${T}] step1 start $(date +%H:%M)"
  nice -n 15 "$SAIGE_RUNNER" "$SAIGE_STEP1_SCRIPT" \
    --phenoFile="$P" --phenoCol="y" \
    --covarColList="age,sex_male,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10" \
    --sampleIDColinphenoFile="sample_id" \
    --traitType="binary" \
    --outputPrefix="$OUT/${T}_v2" \
    --nThreads=12 \
    --plinkFile="$PLINK" \
    --useSparseGRMtoFitNULL=TRUE \
    --sparseGRMFile="$GRM" --sparseGRMSampleIDFile="$GRM_IDS" \
    --skipVarianceRatioEstimation=FALSE \
    --isCateVarianceRatio=TRUE \
    --cateVarRatioMinMACVecExclude="10,20.5" --cateVarRatioMaxMACVecInclude="20.5" \
    --IsOverwriteVarianceRatioFile=TRUE --LOCO=FALSE \
    > "$OUT/${T}_v2.step1.log" 2>&1 || { echo "FAIL: step1 $T"; tail -25 "$OUT/${T}_v2.step1.log"; exit 5; }
  test -s "$OUT/${T}_v2.rda" || { echo "FAIL: $T rda missing"; exit 5; }
  VRLINES=$(wc -l < "$OUT/${T}_v2.varianceRatio.txt")
  [ "$VRLINES" -ge 1 ] || { echo "FAIL: $T variance ratio empty"; exit 5; }
  echo "[${T}] done — VR lines=$VRLINES"
done
grep -H . $OUT/*_v2.varianceRatio.txt | head -20
echo "STAGE3 COMPLETE"
