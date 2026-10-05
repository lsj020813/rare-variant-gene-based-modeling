#!/usr/bin/env bash
: "${CONTAINER_STATE_DIR:?Set CONTAINER_STATE_DIR}"
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
: "${SAIGE_IMAGE:?Set SAIGE_IMAGE}"
: "${UDOCKER_BIN:?Set UDOCKER_BIN}"
export UDOCKER_DIR=${CONTAINER_STATE_DIR}
UD="PROOT_NO_SECCOMP=1 ${UDOCKER_BIN} run --volume=/data:/data ${SAIGE_IMAGE}"
cond_run(){
  local TAG=$1 CONDID=$2
  eval $UD step2_SPAtests.R \
    --vcfFile=${PROJECT_ROOT}/work/run_cond/SORT1_1/mini_id.vcf.gz --vcfFileIndex=${PROJECT_ROOT}/work/run_cond/SORT1_1/mini_id.vcf.gz.csi \
    --vcfField=DS --chrom=1 --AlleleOrder=ref-first --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
    --GMMATmodelFile=${PROJECT_ROOT}/work/ref/saige_step1_v4/tchl_v4.rda \
    --varianceRatioFile=${PROJECT_ROOT}/work/ref/saige_step1_v4/tchl_v4.varianceRatio.txt \
    --groupFile=${PROJECT_ROOT}/work/run_cond/SORT1_1/group.txt --annotation_in_groupTest=all --maxMAF_in_groupTest=0.01 \
    --condition=$CONDID --SAIGEOutputFile=${PROJECT_ROOT}/work/run_cond/SORT1_1/tchl.$TAG > ${PROJECT_ROOT}/work/run_cond/SORT1_1/tchl.$TAG.log 2>&1
  echo \"$TAG done rc=$?\"
}
cond_run diag_lead '1:109817590_G/T'
cond_run diag_ctrl '1:109444458_G/T'
echo DIAG_DONE
