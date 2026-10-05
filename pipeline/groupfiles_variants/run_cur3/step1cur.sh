#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -euo pipefail
OUT=${PROJECT_ROOT}/work/ref/saige_step1
GRM=${GRM_FILE:?Set GRM_FILE}
GRM_IDS=$GRM.sampleIDs.txt
PLINK=${ANALYSIS_ROOT:?Set ANALYSIS_ROOT}/work/saige_gene/varratio_plink_chr22_gp90/chr22.GP90.varratio_mac10plus_retry3
mkdir -p "$OUT"
export UDOCKER_DIR=${CONTAINER_STATE_DIR:?Set CONTAINER_STATE_DIR}
for T in frac1 gastro arth aller per fliv pol oste gasulcer thy cata mi gb bph asth duoulcer; do
  P=${PROJECT_ROOT}/work/ref/pheno_cur/${T}_cur.tsv
  test -s "$P"
  if [ -s "$OUT/${T}_cur.rda" ] && [ -s "$OUT/${T}_cur.varianceRatio.txt" ]; then echo "[${T}] done, skip"; continue; fi
  echo "[${T}] step1 start $(date +%H:%M)"
  PROOT_NO_SECCOMP=1 nice -n 15 ${UDOCKER_BIN:-udocker} run --volume=/data:/data ${SAIGE_IMAGE:?Set SAIGE_IMAGE} \
    step1_fitNULLGLMM.R \
    --phenoFile="$P" --phenoCol="y" \
    --covarColList="age,sex_male,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10" \
    --sampleIDColinphenoFile="sample_id" \
    --traitType="binary" \
    --outputPrefix="$OUT/${T}_cur" \
    --nThreads=4 \
    --plinkFile="$PLINK" \
    --useSparseGRMtoFitNULL=TRUE \
    --sparseGRMFile="$GRM" --sparseGRMSampleIDFile="$GRM_IDS" \
    --skipVarianceRatioEstimation=FALSE \
    --isCateVarianceRatio=TRUE \
    --cateVarRatioMinMACVecExclude="10,20.5" --cateVarRatioMaxMACVecInclude="20.5" \
    --IsOverwriteVarianceRatioFile=TRUE --LOCO=FALSE \
    > "$OUT/${T}_cur.step1.log" 2>&1 || { echo "SKIP: step1 $T did not converge — excluded from curriculum [recorded]"; continue; }
  test -s "$OUT/${T}_cur.rda" || { echo "SKIP: $T rda missing"; continue; }
  VRLINES=$(wc -l < "$OUT/${T}_cur.varianceRatio.txt")
  [ "$VRLINES" -ge 1 ] || { echo "SKIP: $T variance ratio empty"; continue; }
  echo "[${T}] done — VR lines=$VRLINES"
done
grep -H . $OUT/*_cur.varianceRatio.txt | head -20
echo "CURRICULUM_STEP1_COMPLETE"
