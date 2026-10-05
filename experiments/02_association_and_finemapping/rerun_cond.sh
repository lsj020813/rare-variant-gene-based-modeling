#!/usr/bin/env bash
: "${CONTAINER_STATE_DIR:?Set CONTAINER_STATE_DIR}"
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
: "${SAIGE_IMAGE:?Set SAIGE_IMAGE}"
: "${UDOCKER_BIN:?Set UDOCKER_BIN}"
set -uo pipefail
R=${PROJECT_ROOT}/work/ref
W=${PROJECT_ROOT}/work/run_cond
export UDOCKER_DIR=${CONTAINER_STATE_DIR}
run_locus(){
  local L=$1
  local D=$W/$L
  [ -s $D/lead.json ] || { echo "MISS $L lead.json"; return 1; }
  local LEAD=$(python3 -c "import json;l=json.load(open('$D/lead.json'))['lead'];print(l.replace('_',':').replace('/',':'))")
  local CH=${LEAD%%:*}
  PROOT_NO_SECCOMP=1 ${UDOCKER_BIN} run --volume=/data:/data ${SAIGE_IMAGE} \
    step2_SPAtests.R \
    --vcfFile=$D/mini.vcf.gz --vcfFileIndex=$D/mini.vcf.gz.csi \
    --vcfField=DS --chrom=$CH --AlleleOrder=ref-first --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
    --GMMATmodelFile=$R/saige_step1_v4/tchl_v4.rda \
    --varianceRatioFile=$R/saige_step1_v4/tchl_v4.varianceRatio.txt \
    --groupFile=$D/group.txt --annotation_in_groupTest=all --maxMAF_in_groupTest=0.01 \
    --condition=$LEAD \
    --SAIGEOutputFile=$D/tchl.condX > $D/tchl.condX.log 2>&1
  local rc=$?
  if [ ! -s $D/tchl.condX ]; then echo "FAIL $L rc=$rc empty"; return 1; fi
  local ch=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="Pvalue")p=i;if($i=="Pvalue_cond")c=i}} NR>1&&$p!=$c{n++} END{print n+0}' $D/tchl.condX)
  echo "OK $L rc=$rc changed_genes=$ch lead=$LEAD"
}
export -f run_locus; export R W UDOCKER_DIR
printf '%s\n' APOE19 LDLR19 SORT1_1 CETP16 CHR5 APOB2 APO11 | xargs -P 4 -I{} bash -c 'run_locus {}'
echo COND_RERUN_COMPLETE
