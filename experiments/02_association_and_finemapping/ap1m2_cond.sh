#!/usr/bin/env bash
: "${CONTAINER_STATE_DIR:?Set CONTAINER_STATE_DIR}"
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
: "${SAIGE_IMAGE:?Set SAIGE_IMAGE}"
: "${UDOCKER_BIN:?Set UDOCKER_BIN}"
set -uo pipefail
W=${PROJECT_ROOT}/work/run_cond/AP1M2_only; R=${PROJECT_ROOT}/work/ref
mkdir -p $W; cd $W
BCF=bcftools
awk '$1 ~ /^ENSG00000129354/' $R/groupfiles_bwg/chr19.B_3kb_re2g.txt > group.txt
nv=$(awk '$2=="var"{print NF-2}' group.txt)
echo "group vars: $nv"
[ "$nv" -eq 73 ] || { echo "GATE FAIL: vars $nv != 73"; exit 6; }
M=${PROJECT_ROOT}/work/run_cond/LDLR19/mini.vcf.gz
$BCF query -f '%CHROM:%POS:%REF:%ALT\n' $M > mini_keys.txt
miss=0
for k in $(awk '$2=="var"{for(i=3;i<=NF;i++) print $i}' group.txt); do
  grep -qxF "$k" mini_keys.txt || miss=$((miss+1)); done
grep -qF "19:11242307" mini_keys.txt && lead_in=Y || lead_in=N
echo "missing keys in mini: $miss / 73   lead present: $lead_in"
[ "$miss" -eq 0 ] && [ "$lead_in" = Y ] || { echo "GATE FAIL: mini VCF coverage"; exit 6; }
export UDOCKER_DIR=${CONTAINER_STATE_DIR}
run(){
  PROOT_NO_SECCOMP=1 ${UDOCKER_BIN} run --volume=/data:/data ${SAIGE_IMAGE} \
    step2_SPAtests.R --vcfFile=$M --vcfFileIndex=$M.csi \
    --vcfField=DS --chrom=19 --AlleleOrder=ref-first --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
    --GMMATmodelFile=$R/saige_step1_v4/tchl_v4.rda \
    --varianceRatioFile=$R/saige_step1_v4/tchl_v4.varianceRatio.txt \
    --groupFile=$W/group.txt --annotation_in_groupTest=all --maxMAF_in_groupTest=0.01 \
    $2 --SAIGEOutputFile=$W/$1 > $W/$1.log 2>&1
}
run uncond ""
run cond "--condition=19:11242307:G:C"
for f in uncond cond; do
  [ -s $W/$f ] || { echo "GATE FAIL: $f empty"; exit 6; }
  echo "$f rows: $(( $(wc -l < $W/$f) - 1 ))"; done
pu=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Pvalue") p=i; next} {print $p}' $W/uncond)
pc=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Pvalue_cond") p=i; next} {print $p}' $W/cond)
echo "AP1M2 uncond p=$pu   cond p=$pc"
[ "$pu" != "$pc" ] || { echo "GATE FAIL: conditioning had no effect"; exit 6; }
echo "AP1M2_COND_DONE"
