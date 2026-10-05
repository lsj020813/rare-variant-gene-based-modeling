#!/usr/bin/env bash
: "${CONTAINER_STATE_DIR:?Set CONTAINER_STATE_DIR}"
: "${GENOTYPE_DIR:?Set GENOTYPE_DIR}"
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
: "${SAIGE_IMAGE:?Set SAIGE_IMAGE}"
: "${UDOCKER_BIN:?Set UDOCKER_BIN}"
set -euo pipefail
W=${PROJECT_ROOT}/work/run_cond/LDLR2nd; R=${PROJECT_ROOT}/work/ref
M=${PROJECT_ROOT}/work/run_cond/LDLR19/mini.vcf.gz
BCF=bcftools
RAW=${GENOTYPE_DIR}/chr19.vcf.gz
cd $W
L2="19:11343795:A:G"; CH2=19; POS2=11343795
$BCF view -r $CH2:$POS2-$POS2 $RAW -Oz -o lead2.vcf.gz --threads 2
$BCF index -f -c lead2.vcf.gz
n=$($BCF index -n lead2.vcf.gz)
echo "lead2 records: $n"
[ "$n" -ge 1 ] || { echo "GATE FAIL: lead2 extract empty"; exit 6; }
$BCF concat -a $M lead2.vcf.gz -Oz -o mini2.vcf.gz
$BCF sort -T ${PROJECT_ROOT}/work/tmp mini2.vcf.gz -Oz -o mini2.s.vcf.gz
$BCF index -f -c mini2.s.vcf.gz
echo "mini2 records: $($BCF index -n mini2.s.vcf.gz)"
inlead=$($BCF query -r $CH2:$POS2-$POS2 -f '%CHROM:%POS:%REF:%ALT\n' mini2.s.vcf.gz | head -2)
echo "lead2 in mini2: $inlead"
export UDOCKER_DIR=${CONTAINER_STATE_DIR}
for g in ENSG00000079805 ENSG00000142453 ENSG00000127616; do
  d=${PROJECT_ROOT}/work/run_cond/Barm12/$g
  PROOT_NO_SECCOMP=1 ${UDOCKER_BIN} run --volume=/data:/data ${SAIGE_IMAGE} \
    step2_SPAtests.R --vcfFile=$W/mini2.s.vcf.gz --vcfFileIndex=$W/mini2.s.vcf.gz.csi \
    --vcfField=DS --chrom=19 --AlleleOrder=ref-first --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
    --GMMATmodelFile=$R/saige_step1_v4/tchl_v4.rda \
    --varianceRatioFile=$R/saige_step1_v4/tchl_v4.varianceRatio.txt \
    --groupFile=$d/group.txt --annotation_in_groupTest=all --maxMAF_in_groupTest=0.01 \
    --condition=19:11242307:G:C,$L2 \
    --SAIGEOutputFile=$W/$g.cond2c > $W/$g.cond2c.log 2>&1 || true
  [ -s $W/$g.cond2c ] || { echo "FAIL $g: empty"; continue; }
  p1=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Pvalue_cond") p=i; next}{print $p}' $d/cond)
  p2=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Pvalue_cond") p=i; next}{print $p}' $W/$g.cond2c)
  pu=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Pvalue") p=i; next}{print $p}' $d/uncond)
  if [ "$p2" = "$pu" ] || [ "$p2" = "$p1" ]; then echo "INVALID $g (uncond=$pu l1=$p1 l12=$p2)";
  else echo "OK $g  uncond=$pu  lead1=$p1  lead1+2=$p2"; fi
done
echo "LDLR2ND_V3_DONE"
