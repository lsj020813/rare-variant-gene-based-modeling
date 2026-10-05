#!/usr/bin/env bash
: "${CONTAINER_STATE_DIR:?Set CONTAINER_STATE_DIR}"
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
: "${SAIGE_IMAGE:?Set SAIGE_IMAGE}"
: "${UDOCKER_BIN:?Set UDOCKER_BIN}"
W=${PROJECT_ROOT}/work/run_cond/SORT1_1
R=${PROJECT_ROOT}/work/ref
export UDOCKER_DIR=${CONTAINER_STATE_DIR}
PROOT_NO_SECCOMP=1 ${UDOCKER_BIN} run --volume=/data:/data ${SAIGE_IMAGE} \
  step2_SPAtests.R \
  --vcfFile=$W/mini.vcf.gz --vcfFileIndex=$W/mini.vcf.gz.csi \
  --vcfField=DS --chrom=1 --AlleleOrder=ref-first --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  --GMMATmodelFile=$R/saige_step1_v4/tchl_v4.rda \
  --varianceRatioFile=$R/saige_step1_v4/tchl_v4.varianceRatio.txt \
  --groupFile=$W/group.txt --annotation_in_groupTest=all --maxMAF_in_groupTest=0.01 \
  --condition=1:109817590:G:T \
  --SAIGEOutputFile=$W/tchl.colon
echo COLON_DONE rc=$?
