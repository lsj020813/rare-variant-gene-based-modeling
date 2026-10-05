#!/usr/bin/env bash
: "${CONTAINER_STATE_DIR:?Set CONTAINER_STATE_DIR}"
: "${GENOTYPE_DIR:?Set GENOTYPE_DIR}"
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
: "${SAIGE_IMAGE:?Set SAIGE_IMAGE}"
: "${UDOCKER_BIN:?Set UDOCKER_BIN}"
set -uo pipefail
W=${PROJECT_ROOT}/work/run_cond/LDLR2nd; R=${PROJECT_ROOT}/work/ref
M=${PROJECT_ROOT}/work/run_cond/LDLR19/mini.vcf.gz
BCF=bcftools
cd $W
RAW=${GENOTYPE_DIR}/chr19.vcf.gz
$BCF view -r 19:10500000-11500000 -i 'INFO/MAF>0.01 && INFO/R2>0.8' $RAW -Oz -o common_win.vcf.gz --threads 4 2>/dev/null
$BCF index -c common_win.vcf.gz
nc=$($BCF index -n common_win.vcf.gz)
echo "common variants in window: $nc"
[ "$nc" -gt 10 ] || { echo "GATE FAIL: too few common variants"; exit 6; }
export UDOCKER_DIR=${CONTAINER_STATE_DIR}
PROOT_NO_SECCOMP=1 ${UDOCKER_BIN} run --volume=/data:/data ${SAIGE_IMAGE} \
  step2_SPAtests.R --vcfFile=$W/common_win.vcf.gz --vcfFileIndex=$W/common_win.vcf.gz.csi \
  --vcfField=DS --chrom=19 --AlleleOrder=ref-first --minMAF=0.01 --minMAC=20 --LOCO=FALSE \
  --GMMATmodelFile=$R/saige_step1_v4/tchl_v4.rda \
  --varianceRatioFile=$R/saige_step1_v4/tchl_v4.varianceRatio.txt \
  --condition=19:11242307:G:C \
  --SAIGEOutputFile=$W/sv2 > $W/sv2.log 2>&1
[ -s $W/sv2 ] || { echo "GATE FAIL: sv2 empty"; exit 6; }
echo "sv2 rows: $(( $(wc -l < $W/sv2) - 1 ))"
read L2 P2 <<< $(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="CHR")c=i; if($i=="POS")o=i; if($i=="Allele1")a1=i; if($i=="Allele2")a2=i; if($i=="p.value.c"||$i=="p.value_c"||$i=="Pvalue_cond")p=i}; next}
  $o!=11242307 {if(best==""||$p<bp){bp=$p; best=$c":"$o":"$a1":"$a2}} END{print best, bp}' $W/sv2)
echo "lead2: $L2  cond1-p: $P2"
[ -n "$L2" ] && [ "$L2" != ":::" ] || { echo "GATE FAIL: lead2 empty"; exit 6; }
$BCF view -r ${L2%:*:*}:$(echo $L2|cut -d: -f2)-$(echo $L2|cut -d: -f2) $RAW -Oz -o lead2.vcf.gz 2>/dev/null
$BCF index -c lead2.vcf.gz
$BCF concat -a $M lead2.vcf.gz -Oz -o mini2.vcf.gz 2>/dev/null && $BCF sort mini2.vcf.gz -Oz -o mini2.s.vcf.gz 2>/dev/null
$BCF index -c mini2.s.vcf.gz
for g in ENSG00000079805 ENSG00000142453 ENSG00000127616; do
  d=${PROJECT_ROOT}/work/run_cond/Barm12/$g
  PROOT_NO_SECCOMP=1 ${UDOCKER_BIN} run --volume=/data:/data ${SAIGE_IMAGE} \
    step2_SPAtests.R --vcfFile=$W/mini2.s.vcf.gz --vcfFileIndex=$W/mini2.s.vcf.gz.csi \
    --vcfField=DS --chrom=19 --AlleleOrder=ref-first --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
    --GMMATmodelFile=$R/saige_step1_v4/tchl_v4.rda \
    --varianceRatioFile=$R/saige_step1_v4/tchl_v4.varianceRatio.txt \
    --groupFile=$d/group.txt --annotation_in_groupTest=all --maxMAF_in_groupTest=0.01 \
    --condition=19:11242307:G:C,$L2 \
    --SAIGEOutputFile=$W/$g.cond2b > $W/$g.cond2b.log 2>&1
  [ -s $W/$g.cond2b ] || { echo "FAIL $g: empty"; continue; }
  p1=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Pvalue_cond") p=i; next}{print $p}' ${PROJECT_ROOT}/work/run_cond/Barm12/$g/cond)
  p2=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Pvalue_cond") p=i; next}{print $p}' $W/$g.cond2b)
  pu=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Pvalue") p=i; next}{print $p}' ${PROJECT_ROOT}/work/run_cond/Barm12/$g/uncond)
  if [ "$p2" = "$pu" ] || [ "$p2" = "$p1" ]; then echo "INVALID $g: p2 identical (uncond=$pu l1=$p1 l12=$p2)"; else
  echo "OK $g  uncond=$pu  lead1=$p1  lead1+2=$p2"; fi
done
echo "LDLR2ND_V2_DONE"
