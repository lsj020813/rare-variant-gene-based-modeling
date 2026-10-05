#!/usr/bin/env bash
: "${CONTAINER_STATE_DIR:?Set CONTAINER_STATE_DIR}"
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
: "${SAIGE_IMAGE:?Set SAIGE_IMAGE}"
: "${UDOCKER_BIN:?Set UDOCKER_BIN}"
set -uo pipefail
W=${PROJECT_ROOT}/work/run_cond/Barm12; R=${PROJECT_ROOT}/work/ref
mkdir -p $W; cd $W
BCF=bcftools
export UDOCKER_DIR=${CONTAINER_STATE_DIR}
cat > genes.tsv <<'TSV'
ENSG00000143126	1	SORT1_1	1:109817590:G:T
ENSG00000164221	5	CHR5	5:74648603:A:T
ENSG00000152359	5	CHR5	5:74648603:A:T
ENSG00000122008	5	CHR5	5:74648603:A:T
ENSG00000113163	5	CHR5	5:74648603:A:T
ENSG00000130202	19	APOE19	19:45412079:C:T
ENSG00000130204	19	APOE19	19:45412079:C:T
ENSG00000104972	19	APOE19	19:45412079:C:T
ENSG00000069399	19	APOE19	19:45412079:C:T
ENSG00000079805	19	LDLR19	19:11242307:G:C
ENSG00000142453	19	LDLR19	19:11242307:G:C
ENSG00000127616	19	LDLR19	19:11242307:G:C
TSV
run_gene(){
  g=$1; ch=$2; loc=$3; lead=$4
  d=$W/$g; mkdir -p $d
  awk -v G=$g '$1 ~ "^"G' $R/groupfiles_bwg/chr$ch.B_3kb_re2g.txt > $d/group.txt
  nv=$(awk '$2=="var"{print NF-2}' $d/group.txt)
  M=${PROJECT_ROOT}/work/run_cond/$loc/mini.vcf.gz
  $BCF query -f '%CHROM:%POS:%REF:%ALT\n' $M > $d/mini_keys.txt
  miss=0
  for k in $(awk '$2=="var"{for(i=3;i<=NF;i++) print $i}' $d/group.txt); do
    grep -qxF "$k" $d/mini_keys.txt || miss=$((miss+1)); done
  if [ "$miss" -gt 0 ]; then echo "SKIP $g: $miss/$nv keys not in mini"; return; fi
  for mode in uncond cond; do
    extra=""; [ $mode = cond ] && extra="--condition=$lead"
    PROOT_NO_SECCOMP=1 ${UDOCKER_BIN} run --volume=/data:/data ${SAIGE_IMAGE} \
      step2_SPAtests.R --vcfFile=$M --vcfFileIndex=$M.csi \
      --vcfField=DS --chrom=$ch --AlleleOrder=ref-first --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
      --GMMATmodelFile=$R/saige_step1_v4/tchl_v4.rda \
      --varianceRatioFile=$R/saige_step1_v4/tchl_v4.varianceRatio.txt \
      --groupFile=$d/group.txt --annotation_in_groupTest=all --maxMAF_in_groupTest=0.01 \
      $extra --SAIGEOutputFile=$d/$mode > $d/$mode.log 2>&1
    [ -s $d/$mode ] || { echo "FAIL $g $mode: empty"; return; }
  done
  pu=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Pvalue") p=i; next}{print $p}' $d/uncond)
  pc=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Pvalue_cond") p=i; next}{print $p}' $d/cond)
  [ "$pu" != "$pc" ] || { echo "FAIL $g: conditioning no effect"; return; }
  echo "OK $g nv=$nv uncond=$pu cond=$pc"
}
export -f run_gene; export W R BCF UDOCKER_DIR
cat genes.tsv | xargs -P 4 -n 4 bash -c 'run_gene "$0" "$1" "$2" "$3"' >> $W/results.log 2>&1
echo "BARM12_DONE" >> $W/results.log
