#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT to the working data root}"
[[ "$PROJECT_ROOT" =~ [^[:space:]] ]] || { printf "%s\n" "PROJECT_ROOT must be nonblank" >&2; exit 2; }
CODE_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
export CODE_ROOT
set -uo pipefail
cd ${PROJECT_ROOT}/work/prs

startwd() {
  tmux kill-session -t wd_$3 2>/dev/null
  tmux new-session -d -s wd_$3 "bash \"${CODE_ROOT}/infra/$1\" '$2' $3"; sleep 5
  pgrep -f "$1 .*$3" >/dev/null && [ -s logs/wd_$3.log ] || { echo "WATCHDOG FAILED $3 -> abort"; exit 9; }
  echo "watchdog $3 ok $(date -Is)"
}

while [ ! -s out/hei/chain_linear.done ]; do sleep 120; done
R=${PROJECT_ROOT}/work/ref; OUT=$R/saige_step1_hei; mkdir -p $OUT
GRM="${SPARSE_GRM_PATH:?Set SPARSE_GRM_PATH}"
PLINK=$R/vr_plink_v3/vr_v3
O=$OUT/hei_v1; [ -s $O.rda ] && [ -s $O.varianceRatio.txt ] && { echo skip; exit 0; }
startwd wd_prs.sh '[s]tep1_fitNULLGLMM.R' heis1b
PROOT_NO_SECCOMP=1 "${UDOCKER_BIN:-udocker}" run --volume=/data:/data "${SAIGE_IMAGE:?Set SAIGE_IMAGE}" step1_fitNULLGLMM.R \
  --sparseGRMFile="$GRM" --sparseGRMSampleIDFile="$GRM.sampleIDs.txt" --useSparseGRMtoFitNULL=TRUE --plinkFile="$PLINK" \
  --phenoFile="$R/pheno_hei/hei_rint.tsv" --phenoCol="y" --sampleIDColinphenoFile="sample_id" \
  --covarColList="${PHENO_COVARIATE_COLUMNS:?Set PHENO_COVARIATE_COLUMNS}" --traitType=quantitative --invNormalize=FALSE \
  --isCateVarianceRatio=TRUE --cateVarRatioMinMACVecExclude="10,20.5" --cateVarRatioMaxMACVecInclude="20.5" \
  --skipVarianceRatioEstimation=FALSE --IsOverwriteVarianceRatioFile=TRUE --LOCO=FALSE --nThreads=12 --outputPrefix="$O" > $O.step1.log 2>&1
rc=$?; [ $rc -eq 0 ] && [ -s $O.rda ] && [ -s $O.varianceRatio.txt ] && echo "ok $(date -Is)" > $OUT/hei_v1.done; tmux kill-session -t wd_heis1b; echo "step1 rc=$rc"
