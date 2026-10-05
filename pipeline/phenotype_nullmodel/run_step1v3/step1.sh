#!/usr/bin/env bash
set -uo pipefail
: "${PROJECT_ROOT:?Set PROJECT_ROOT to a generic workspace root}"
RUN_DIR=${RUN_DIR:-"${PROJECT_ROOT}/work/run_step1v3"}
mkdir -p "$RUN_DIR"
exec 9>"$RUN_DIR/.lock"
flock -n 9 || { echo "another instance running"; exit 0; }

R=${REFERENCE_DIR:-"${PROJECT_ROOT}/work/ref"}
P=${PHENO_V3_DIR:-"$R/pheno_v3"}
OUT=$R/saige_step1_v3
GRM=${GRM_FILE:?Set GRM_FILE}
GRM_IDS=${GRM_SAMPLE_IDS_FILE:?Set GRM_SAMPLE_IDS_FILE}
PLINK=${PLINK_PREFIX:?Set PLINK_PREFIX}
COVAR="age,sex_male,CT,NC,PC1,PC2,PC3,PC4,PC5"
mkdir -p "$OUT"
NCORE=$(nproc)
NTHREADS=${NTHREADS:-$(( NCORE / 8 ))}
[ "$NTHREADS" -lt 4 ] && NTHREADS=4
echo "cores=$NCORE nThreads=$NTHREADS"
SAIGE_RUNNER=${SAIGE_RUNNER:-Rscript}
SAIGE_STEP1_SCRIPT=${SAIGE_STEP1_SCRIPT:?Set SAIGE_STEP1_SCRIPT}

for f in "$GRM" "$GRM_IDS" "$PLINK.bed" "$PLINK.bim" "$PLINK.fam"; do
  [ -s "$f" ] || { echo "PREFLIGHT FAIL missing $f"; exit 3; }
done
for T in htn dm lip tchl; do
  [ -s "$P/${T}_v3.tsv" ] || { echo "PREFLIGHT FAIL missing $P/${T}_v3.tsv"; exit 3; }
done
echo "preflight OK"

fit() {
  local T=$1 TYPE=$2
  local O="$OUT/${T}_v3"
  local INVNORM=""
  [ "$TYPE" = "quantitative" ] && INVNORM="--invNormalize=FALSE"
  [ -s "$O.rda" ] && [ -s "$O.varianceRatio.txt" ] && { echo "skip $T"; return 0; }
  "$SAIGE_RUNNER" "$SAIGE_STEP1_SCRIPT" \
    --sparseGRMFile="$GRM" --sparseGRMSampleIDFile="$GRM_IDS" \
    --useSparseGRMtoFitNULL=TRUE \
    --plinkFile="$PLINK" \
    --phenoFile="$P/${T}_v3.tsv" --phenoCol="y" --sampleIDColinphenoFile="sample_id" \
    --covarColList="$COVAR" \
    --traitType="$TYPE" \
    $INVNORM \
    --isCateVarianceRatio=TRUE \
    --cateVarRatioMinMACVecExclude="10,20.5" --cateVarRatioMaxMACVecInclude="20.5" \
    --skipVarianceRatioEstimation=FALSE --IsOverwriteVarianceRatioFile=TRUE \
    --LOCO=FALSE --nThreads="$NTHREADS" \
    --outputPrefix="$O" > "$O.step1.log" 2>&1
  local rc=$?
  if [ "$rc" -ne 0 ] || [ ! -s "$O.rda" ] || [ ! -s "$O.varianceRatio.txt" ]; then
    echo "FAIL $T rc=$rc rda=$( [ -s "$O.rda" ] && echo Y || echo N ) vr=$( [ -s "$O.varianceRatio.txt" ] && echo Y || echo N )"
    tail -5 "$O.step1.log"
    return 1
  fi
  echo "OK $T  vr=$(head -1 "$O.varianceRatio.txt" | awk '{print $1}')"
}

fit htn  binary
fit dm   binary
fit lip  binary
fit tchl quantitative

n=$(ls "$OUT"/*_v3.rda 2>/dev/null | wc -l)
echo "null models produced: $n/4"
[ "$n" -eq 4 ] || { echo "GATE FAIL: expected 4 null models"; exit 6; }
echo "STEP1_V3_COMPLETE"
