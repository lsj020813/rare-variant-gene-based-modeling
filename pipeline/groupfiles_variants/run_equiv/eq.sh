#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
D=${PROJECT_ROOT}/work/run_equiv
GF=${PROJECT_ROOT}/work/ref/groupfiles_v1/chr21.groupfile_V1.txt
BAND=${PROJECT_ROOT}/work/ref/band_vcf/chr21.band.vcf.gz
PILOT=${ANALYSIS_ROOT:?Set ANALYSIS_ROOT}
MODEL=$PILOT/work/saige_gene/grch38_chr22/step1/tchl_primary.rda
VR=$PILOT/work/saige_gene/grch38_chr22/step1/tchl_primary_with_vr_chr22_gp90_retry3.varianceRatio.txt
mkdir -p $D; cd $D
export UDOCKER_DIR=${CONTAINER_STATE_DIR:?Set CONTAINER_STATE_DIR}

awk '{ if ($1 != prev) { n++; prev=$1 } if (n<=60) print }' $GF > A_60.txt
awk '{ if ($1 != prev) { n++; prev=$1 } if (n<=30) print }' $GF > B_first30.txt
awk '{ if ($1 != prev) { n++; prev=$1 } if (n>30 && n<=60) print }' $GF > B_second30.txt
for f in A_60 B_first30 B_second30; do
  printf '%s: %s genes, %s lines\n' $f $(awk '{print $1}' $f.txt | sort -u | wc -l) $(wc -l < $f.txt)
done
cat B_first30.txt B_second30.txt | sort > /tmp/b.s; sort A_60.txt > /tmp/a.s
cmp -s /tmp/a.s /tmp/b.s || { echo "GATE FAIL: B union != A"; exit 5; }
echo "GATE OK: B union == A byte-identical"

run(){
  local TAG=$1 G=$2; shift 2
  PROOT_NO_SECCOMP=1 nice -n 5 ${UDOCKER_BIN:-udocker} run --volume=/data:/data ${SAIGE_IMAGE:?Set SAIGE_IMAGE} \
    step2_SPAtests.R \
    --vcfFile="$BAND" --vcfFileIndex="$BAND.csi" --vcfField="DS" --chrom=21 \
    --minMAF=0 --minMAC=0.5 --maxMAF_in_groupTest=0.01 \
    --GMMATmodelFile="$MODEL" --varianceRatioFile="$VR" \
    --SAIGEOutputFile="$D/$TAG" --groupFile="$G" \
    --annotation_in_groupTest="all" --is_fastTest=TRUE \
    --LOCO=FALSE --AlleleOrder=ref-first "$@" > $D/$TAG.log 2>&1 \
    && echo "[$TAG] rows=$(wc -l < $D/$TAG)" || { echo "[$TAG] FAIL"; tail -3 $D/$TAG.log; }
}
run A     $D/A_60.txt          &
run B1    $D/B_first30.txt     &
run B2    $D/B_second30.txt    &
run C     $D/A_60.txt --groups_per_chunk=7 &
wait
echo EQUIV_RUNS_DONE
