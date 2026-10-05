#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -o pipefail
T=$1; RID=$2; ARMS=${3:-gate,phi,perm}; NPERM=${4:-20}; TAG=${5:-}
W=${PROJECT_ROOT}/work/run_trackA; F=${PROJECT_ROOT}/work/run_ourfm/fm/$T
export TMPDIR=${PROJECT_ROOT}/work/tmp TMP=${PROJECT_ROOT}/work/tmp TEMP=${PROJECT_ROOT}/work/tmp
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
export FM_MAXVAR_MEM=${FM_MAXVAR_MEM:-45000}
N=$(echo $RID | sed 's/.*_chr\([0-9]*\)_.*/\1/')
D=$W/fm/$RID${TAG:+.$TAG}; mkdir -p $W/fm
[ -s $W/fm/$RID.gate.done ] && ARMS=$(echo "$ARMS" | sed "s/gate,//; s/,gate//; s/^gate$//")
[ -n "$ARMS" ] || { echo "[$RID] arms 비어있음"; exit 0; }
[ -s $D.done ] && { echo "[$RID] 이미 완료"; exit 0; }
PHI=$W/phi/$RID.phi.tsv; [ "$ARMS" = "gate" ] && PHI=-
for f in $F/$RID.step2.txt.done $F/$RID.ld.done $F/$RID.ld.bin $F/$RID.ld.vars $F/$RID.pip.tsv; do [ -s $f ] || { echo "[$RID] GATE FAIL: 입력 없음 $f"; exit 2; }; done
NV=$(wc -l < $F/$RID.ld.vars); SZ=$(stat -c %s $F/$RID.ld.bin); [ "$SZ" -eq $(( NV * NV * 4 )) ] || { echo "[$RID] GATE FAIL: LD bin 크기"; exit 3; }
if [ "$PHI" != "-" ]; then [ -s $PHI ] || { echo "[$RID] GATE FAIL: φ 없음"; exit 2; }; NP=$(( $(wc -l < $PHI) - 1 )); [ "$NP" -eq "$NV" ] || { echo "[$RID] GATE FAIL: φ 행 $NP != 변이 $NV"; exit 3; }; fi
T0=$(date +%s)
FM_UNIFORM_PIP=$F/$RID.pip.tsv nice -n 10 Rscript $W/scripts/fm_susie_prior2.R $T $N $RID $F/$RID.step2.txt $F/$RID.ld $D $PHI $ARMS $NPERM > $D.prior.log 2>&1
RC=$?
if grep -q '^PRIOR_SKIP' $D.prior.log; then echo "[$RID] SKIP $(grep '^PRIOR_SKIP' $D.prior.log)"; grep '^PRIOR_SKIP' $D.prior.log > $D.skip; exit 0; fi
grep -q '^PRIOR_OK' $D.prior.log && [ -s $D.arms.tsv ] || { echo "[$RID] GATE FAIL: rc=$RC"; tail -3 $D.prior.log; exit 8; }
echo ok > $D.done
echo "[$RID] $(grep '^GATE' $D.prior.log) | $(grep '^PHI' $D.prior.log) | $(grep '^PRIOR_OK' $D.prior.log) $(( $(date +%s)-T0 ))s"
