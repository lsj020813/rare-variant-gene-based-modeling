#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
W=${PROJECT_ROOT}/work/run_ourfm; cd $W
R=/usr/bin/Rscript
[ -x $R ] || R=Rscript
echo "[chk] $(date '+%m-%d %H:%M') START"
T=tchl; RID=tchl_chr3_r1; D=$W/fm/$T/$RID
if [ -s $D.step2.txt ] && [ -s $D.ld.bin ]; then
  OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 $R $W/scripts/fm_equiv3.R $D.step2.txt $D.ld > $W/logs/equiv3.log 2>&1
  echo "[chk] equiv rc=$? :: $(grep -E 'EQ_RESULT|EQ_PASS|EQ_DIFF' $W/logs/equiv3.log | tr '\n' ' ')"
else
  echo "[chk] equiv SKIP: 캐시 없음 ($D)"
fi
for RID in tchl_chr1_r1 tchl_chr1_r2 tchl_chr1_r3 lip_chr1_r1; do
  T=${RID%%_*}; N=1; D=$W/fm/$T/$RID
  [ -s $D.step2.txt ] && [ -s $D.ld.bin ] || { echo "[chk] $RID SKIP 캐시없음"; continue; }
  OUT=$D.mi1000
  [ -s $OUT.summary.tsv ] && { echo "[chk] $RID 이미 완료"; continue; }
  S=$(date +%s)
  FM_MAX_ITER=1000 FM_MAXVAR_MEM=30000 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 \
    $R $W/scripts/fm_susie.R $T $N $RID $D.step2.txt $D.ld $OUT > $W/logs/mi1000.$RID.log 2>&1
  echo "[chk] $RID rc=$? $(( $(date +%s)-S ))s :: $(grep -E 'SUSIE_OK|GATE FAIL|SUSIE_SKIP|Error' $W/logs/mi1000.$RID.log | head -1 | cut -c1-140)"
done
echo "[chk] $(date '+%m-%d %H:%M') CHK_DONE"
