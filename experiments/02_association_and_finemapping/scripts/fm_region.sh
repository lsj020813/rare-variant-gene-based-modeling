#!/usr/bin/env bash
: "${CONDA_PREFIX:?Set CONDA_PREFIX}"
: "${GENOTYPE_DIR:?Set GENOTYPE_DIR}"
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
: "${SAIGE_IMAGE:?Set SAIGE_IMAGE}"
: "${UDOCKER_BIN:?Set UDOCKER_BIN}"
set -uo pipefail
T=$1; N=$2; S=$3; E=$4; RID=$5
W=${PROJECT_ROOT}/work/run_ourfm; R=${PROJECT_ROOT}/work/ref
export PATH=${CONDA_PREFIX}/bin:$PATH
export TMPDIR=${PROJECT_ROOT}/work/tmp TMP=${PROJECT_ROOT}/work/tmp TEMP=${PROJECT_ROOT}/work/tmp
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
export FM_MAXVAR_MEM=${FM_MAXVAR_MEM:-45000}
ORIG=${GENOTYPE_DIR}/chr$N.vcf.gz
UD=${UDOCKER_BIN}
PLINK=plink2
MAXVAR=${MAXVAR:-45000}
mkdir -p $W/regcache $W/fm/$T
REG=$W/regcache/chr$N.$S.$E; V=$REG.r2ge07.vcf.gz; D=$W/fm/$T/$RID
[ -s $D.done ] && { echo "[$RID] 이미 완료"; exit 0; }
[ -s $D.skip ] && { echo "[$RID] 이미 skip: $(cat $D.skip)"; exit 0; }
[ -s $W/keep/$T.keep ] || { echo "[$RID] GATE FAIL: keep 없음"; exit 2; }
if [ ! -s $V.done ]; then
  T0=$(date +%s)
  bcftools view -r $N:$S-$E -i 'INFO/R2>=0.7' --threads 2 -Ou $ORIG | bcftools annotate -x FORMAT/GT,FORMAT/GP --threads 4 -Oz -o $V.tmp
  RC=$?; [ $RC -eq 0 ] || { echo "[$RID] GATE FAIL: 추출 파이프 rc=$RC"; rm -f $V.tmp; exit 3; }
  EOF_=$(tail -c 28 $V.tmp | od -An -tx1 | tr -d ' \n')
  [ "$EOF_" = "1f8b08040000000000ff0600424302001b0003000000000000000000" ] || { echo "[$RID] GATE FAIL: BGZF EOF 없음"; rm -f $V.tmp; exit 3; }
  mv $V.tmp $V && rm -f $V.tmp.csi
  bcftools index -c -f $V || { echo "[$RID] GATE FAIL: 색인"; exit 3; }
  NREC=$(bcftools index -n $V); [ "${NREC:-0}" -gt 0 ] 2>/dev/null || { echo "[$RID] GATE FAIL: 추출 0건 (index -n='$NREC')"; exit 3; }
  echo "ok $NREC" > $V.done
  echo "[$RID] EXTRACT_DONE n=$NREC $(( $(date +%s)-T0 ))s"
fi
NREC=$(cut -d' ' -f2 $V.done)
[ "$NREC" -le $MAXVAR ] || { echo "[$RID] SKIP_TOO_LARGE n=$NREC > $MAXVAR"; echo "too_large $NREC" > $D.skip; echo "$RID $NREC" >> ${PROJECT_ROOT}/work/run_ourfm/oversize.txt; exit 0; }
ST=$D.step2.txt
if [ ! -s $ST.done ]; then
  T0=$(date +%s)
  $UD run --volume=/data:/data ${SAIGE_IMAGE} step2_SPAtests.R \
    --vcfFile=$V --vcfFileIndex=$V.csi --vcfField=DS --chrom=$N --AlleleOrder=ref-first --minMAF=0 --minMAC=20 \
    --GMMATmodelFile=$R/saige_step1_v4/${T}_v4.rda --varianceRatioFile=$R/saige_step1_v4/${T}_v4.varianceRatio.txt \
    --LOCO=FALSE --SAIGEOutputFile=$ST > $ST.log 2>&1
  RC=$?
  [ -s $ST ] || { echo "[$RID] GATE FAIL: step2 출력 없음 rc=$RC"; tail -3 $ST.log; exit 4; }
  HDR=$(head -1 $ST); case "$HDR" in CHR*POS*MarkerID*) ;; *) echo "[$RID] GATE FAIL: step2 헤더"; exit 5;; esac
  NR_=$(( $(wc -l < $ST) - 1 )); [ $NR_ -gt 1 ] || { echo "[$RID] GATE FAIL: step2 행 $NR_"; exit 6; }
  echo "ok $NR_" > $ST.done; echo "[$RID] STEP2_DONE n=$NR_ of $NREC $(( $(date +%s)-T0 ))s"
fi
LD=$D.ld
if [ ! -s $LD.done ]; then
  T0=$(date +%s)
  LD_THREADS=4 python3 $W/scripts/fm_ld.py $V $W/keep/$T.keep $D $ST > $LD.py.log 2>&1
  RC=$?; grep -q '^LD_OK' $LD.py.log && [ -s $LD.bin ] && [ -s $LD.vars ] || { echo "[$RID] GATE FAIL: LD rc=$RC"; tail -3 $LD.py.log; exit 7; }
  NV=$(wc -l < $LD.vars); SZ=$(stat -c %s $LD.bin)
  [ "$SZ" -eq $(( NV * NV * 4 )) ] || { echo "[$RID] GATE FAIL: LD bin 크기 $SZ != $NV^2*4"; exit 7; }
  echo "ok $NV" > $LD.done; echo "[$RID] LD_DONE n=$NV $(( $(date +%s)-T0 ))s"
fi
T0=$(date +%s)
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 Rscript $W/scripts/fm_susie.R $T $N $RID $ST $LD $D > $D.susie.log 2>&1
RC=$?
if grep -q '^SUSIE_SKIP' $D.susie.log; then echo "[$RID] SKIP_MEM $(grep '^SUSIE_SKIP' $D.susie.log)"; grep '^SUSIE_SKIP' $D.susie.log | sed 's/SUSIE_SKIP //' > $D.skip; exit 0; fi
grep -q '^SUSIE_OK' $D.susie.log && [ -s $D.pip.tsv ] && [ -s $D.summary.tsv ] || { echo "[$RID] GATE FAIL: susie rc=$RC"; tail -3 $D.susie.log; exit 8; }
echo ok > $D.done
echo "[$RID] $(grep '^SUSIE_OK' $D.susie.log) $(( $(date +%s)-T0 ))s"
