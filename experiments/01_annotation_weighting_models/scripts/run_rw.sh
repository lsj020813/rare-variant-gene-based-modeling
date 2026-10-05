#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -o pipefail
T=$1; RID=$2; NPERM=${3:-1000}
W=${PROJECT_ROOT}/work/run_trackA; F=${PROJECT_ROOT}/work/run_ourfm/fm/$T
export TMPDIR=${PROJECT_ROOT}/work/tmp TMP=${PROJECT_ROOT}/work/tmp TEMP=${PROJECT_ROOT}/work/tmp
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 FM_MAXVAR_MEM=${FM_MAXVAR_MEM:-45000}
N=$(echo $RID | sed 's/.*_chr\([0-9]*\)_.*/\1/'); D=$W/rw/$RID; mkdir -p $W/rw
[ -s $D.done ] && { echo "[$RID] 이미 완료"; exit 0; }
for f in $F/$RID.step2.txt.done $F/$RID.ld.done $F/$RID.ld.bin $F/$RID.ld.vars $F/$RID.pip.tsv $W/phi/$RID.phi.tsv; do [ -s $f ] || { echo "[$RID] GATE FAIL: 입력 없음 $f"; exit 2; }; done
NV=$(wc -l < $F/$RID.ld.vars); SZ=$(stat -c %s $F/$RID.ld.bin); [ "$SZ" -eq $(( NV * NV * 4 )) ] || { echo "[$RID] GATE FAIL: LD bin 크기"; exit 3; }
T0=$(date +%s)
FM_UNIFORM_PIP=$F/$RID.pip.tsv nice -n 10 Rscript $W/scripts/fm_reweight.R $T $N $RID $F/$RID.step2.txt $F/$RID.ld $D $W/phi/$RID.phi.tsv $NPERM > $D.rw.log 2>&1
RC=$?
if grep -q '^RW_SKIP' $D.rw.log; then grep '^RW_SKIP' $D.rw.log > $D.skip; echo "[$RID] SKIP"; exit 0; fi
grep -q '^RW_OK' $D.rw.log && [ -s $D.rw.arms.tsv ] && [ -s $D.rw.perms.tsv ] || { echo "[$RID] GATE FAIL: rc=$RC"; tail -3 $D.rw.log; exit 8; }
echo ok > $D.done
echo "[$RID] $(grep '^GATE' $D.rw.log) | $(grep '^PHI_RW' $D.rw.log) | $(grep '^RW_OK' $D.rw.log) $(( $(date +%s)-T0 ))s"
