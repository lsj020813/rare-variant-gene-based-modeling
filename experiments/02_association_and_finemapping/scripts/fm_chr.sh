#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
N=$1; SM=${2:-}
W=${PROJECT_ROOT}/work/run_ourfm; G=${PROJECT_ROOT}/work/ref/gwas05
export TMPDIR=${PROJECT_ROOT}/work/tmp TMP=${PROJECT_ROOT}/work/tmp TEMP=${PROJECT_ROOT}/work/tmp
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
LOG=$W/logs/chr$N.log; exec >> $LOG 2>&1
echo "[chr$N] START $(date '+%F %T') smoke='$SM'"
for i in $(seq 1 1680); do ok=1; for T in tchl htn dm lip; do [ -s $G/$T.chr$N.txt.done ] || ok=0; done; [ $ok -eq 1 ] && break; sleep 30; done
for T in tchl htn dm lip; do [ -s $G/$T.chr$N.txt.done ] || { echo "[chr$N] GATE FAIL: GWAS 마커 대기 초과 ($T)"; exit 2; }; done
echo "[chr$N] GWAS ready $(date '+%F %T')"
AV=$(free -g | awk 'NR==2{print $7}'); [ "$AV" -ge 100 ] || { echo "[chr$N] GATE FAIL: 가용 램 ${AV}G <100G"; exit 9; }
FAIL=0
for T in tchl htn dm lip; do
  python3 $W/scripts/fm_regions.py $T $N $SM || { echo "[chr$N $T] GATE FAIL: regions"; FAIL=1; continue; }
  tail -n +2 $W/regions/$T.chr$N.regions.tsv | while IFS=$'\t' read -r RID CH S E NL MP SMK; do
    AV=$(free -g | awk 'NR==2{print $7}'); [ "$AV" -ge 100 ] || { echo "[$RID] GATE FAIL: 가용 램 ${AV}G <100G — 중단"; exit 9; }
    WA=$(vmstat 1 2 | tail -1 | awk '{print $16}'); [ "${WA:-0}" -le 20 ] || { echo "[$RID] GATE FAIL: iowait ${WA}% >20 — 중단"; exit 9; }
    bash $W/scripts/fm_region.sh $T $CH $S $E $RID; RC=$?; [ $RC -eq 0 ] || { echo "[$RID] REGION FAIL rc=$RC"; exit 1; }
  done; RCL=$?
  [ $RCL -eq 9 ] && { echo "[chr$N] 서버 안전 게이트 — 염색체 중단"; exit 9; }
  [ $RCL -eq 0 ] || { FAIL=1; echo "[chr$N $T] 구역 루프 실패 rc=$RCL"; }
  OUTP=$W/our_pip/$T.chr$N.tsv
  ls $W/fm/$T/${T}_chr${N}_r*.pip.tsv >/dev/null 2>&1 && { awk 'FNR==1 && NR!=1{next} {print}' $W/fm/$T/${T}_chr${N}_r*.pip.tsv > $OUTP; } || : > $OUTP
  NREG=$(( $(wc -l < $W/regions/$T.chr$N.regions.tsv) - 1 )); NDONE=$(ls $W/fm/$T/ 2>/dev/null | grep -cE "^${T}_chr${N}_r[0-9]+\.done$"); NSKIP=$(ls $W/fm/$T/ 2>/dev/null | grep -cE "^${T}_chr${N}_r[0-9]+\.skip$")
  echo "[chr$N $T] TRAIT_SUMMARY regions=$NREG done=$NDONE skip=$NSKIP"
  [ $(( NDONE + NSKIP )) -eq $NREG ] && echo "ok $NDONE $NSKIP" > $OUTP.done
done
[ $FAIL -eq 0 ] && echo "ok" > $W/logs/chr$N.done
echo "[chr$N] CHR_END fail=$FAIL $(date '+%F %T')"
