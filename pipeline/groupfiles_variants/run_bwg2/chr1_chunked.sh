#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
R=${PROJECT_ROOT}/work/ref
O=$R/saige_step2_bwg_chr1
G=$R/groupfiles_bwg_chunks
mkdir -p $O
export UDOCKER_DIR=${CONTAINER_STATE_DIR:?Set CONTAINER_STATE_DIR}
run_one(){
  T=$1; P=$2
  out=$O/$T.chr1.$P
  [ -s $out.done ] && return 0
  PROOT_NO_SECCOMP=1 ${UDOCKER_BIN:-udocker} run --volume=/data:/data ${SAIGE_IMAGE:?Set SAIGE_IMAGE} \
    step2_SPAtests.R \
    --vcfFile=$R/band_vcf/chr1.band.vcf.gz \
    --vcfFileIndex=$R/band_vcf/chr1.band.vcf.gz.csi \
    --vcfField=DS --chrom=1 --AlleleOrder=ref-first \
    --minMAF=0 --minMAC=0.5 --LOCO=FALSE --is_fastTest=TRUE \
    --GMMATmodelFile=$R/saige_step1_v4/${T}_v4.rda \
    --varianceRatioFile=$R/saige_step1_v4/${T}_v4.varianceRatio.txt \
    --groupFile=$G/chr1.$P.txt \
    --annotation_in_groupTest=all --maxMAF_in_groupTest=0.01 \
    --SAIGEOutputFile=$out > $out.log 2>&1
  rc=$?
  exp=$(awk '$2=="var"{print $1}' $G/chr1.$P.txt | sort -u | wc -l)
  nrow=$([ -s $out ] && echo $(( $(wc -l < $out) - 1 )) || echo 0)
  if [ $rc -eq 0 ] && [ "$nrow" -eq "$exp" ]; then echo ok > $out.done; echo "  OK $T $P $nrow/$exp"
  else echo "  FAIL $T $P rc=$rc rows=$nrow/$exp"; fi
}
export -f run_one; export O G R UDOCKER_DIR
{ for T in tchl htn dm lip; do for f in $G/chr1.part*.txt; do
    echo "$T $(basename $f .txt | sed 's/chr1\.//')"; done; done; } | \
  xargs -P 75 -n2 bash -c 'run_one "$@"' _
echo "done markers: $(ls $O/*.done 2>/dev/null | wc -l)/584"
echo "CHR1_CHUNKED_COMPLETE"
