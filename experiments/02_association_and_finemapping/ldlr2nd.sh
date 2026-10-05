#!/usr/bin/env bash
: "${CONTAINER_STATE_DIR:?Set CONTAINER_STATE_DIR}"
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
: "${SAIGE_IMAGE:?Set SAIGE_IMAGE}"
: "${UDOCKER_BIN:?Set UDOCKER_BIN}"
set -uo pipefail
W=${PROJECT_ROOT}/work/run_cond/LDLR2nd; R=${PROJECT_ROOT}/work/ref
M=${PROJECT_ROOT}/work/run_cond/LDLR19/mini.vcf.gz
mkdir -p $W; cd $W
export UDOCKER_DIR=${CONTAINER_STATE_DIR}
PROOT_NO_SECCOMP=1 ${UDOCKER_BIN} run --volume=/data:/data ${SAIGE_IMAGE} \
  step2_SPAtests.R --vcfFile=$M --vcfFileIndex=$M.csi \
  --vcfField=DS --chrom=19 --AlleleOrder=ref-first --minMAF=0.01 --minMAC=20 --LOCO=FALSE \
  --GMMATmodelFile=$R/saige_step1_v4/tchl_v4.rda \
  --varianceRatioFile=$R/saige_step1_v4/tchl_v4.varianceRatio.txt \
  --condition=19:11242307:G:C \
  --SAIGEOutputFile=$W/sv_cond1 > $W/sv_cond1.log 2>&1
[ -s $W/sv_cond1 ] || { echo "FAIL: sv scan empty"; exit 6; }
echo "sv rows: $(( $(wc -l < $W/sv_cond1) - 1 ))"
lead2=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="MarkerID")m=i; if($i=="p.value_c"||$i=="p.value.c"||$i=="Pvalue_cond")p=i; if($i=="AF_Allele2")af=i}; next}
  $m!="19:11242307:G:C" && $af>0.01 && $af<0.99 {if(best==""||$p<bp){bp=$p; best=$m}} END{print best"\t"bp}' $W/sv_cond1)
echo "lead2: $lead2"
L2=$(echo "$lead2" | cut -f1)
[ -n "$L2" ] || { echo "FAIL: no lead2 found"; exit 6; }
for g in ENSG00000079805 ENSG00000142453 ENSG00000127616; do
  d=${PROJECT_ROOT}/work/run_cond/Barm12/$g
  PROOT_NO_SECCOMP=1 ${UDOCKER_BIN} run --volume=/data:/data ${SAIGE_IMAGE} \
    step2_SPAtests.R --vcfFile=$M --vcfFileIndex=$M.csi \
    --vcfField=DS --chrom=19 --AlleleOrder=ref-first --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
    --GMMATmodelFile=$R/saige_step1_v4/tchl_v4.rda \
    --varianceRatioFile=$R/saige_step1_v4/tchl_v4.varianceRatio.txt \
    --groupFile=$d/group.txt --annotation_in_groupTest=all --maxMAF_in_groupTest=0.01 \
    --condition=19:11242307:G:C,$L2 \
    --SAIGEOutputFile=$W/$g.cond2 > $W/$g.cond2.log 2>&1
  [ -s $W/$g.cond2 ] || { echo "FAIL $g cond2 empty"; continue; }
  p1=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Pvalue_cond") p=i; next}{print $p}' ${PROJECT_ROOT}/work/run_cond/Barm12/$g/cond)
  p2=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Pvalue_cond") p=i; next}{print $p}' $W/$g.cond2)
  echo "OK $g  lead1-cond=$p1  lead1+2-cond=$p2"
done
echo "LDLR2ND_DONE"
