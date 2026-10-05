#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
W=${PROJECT_ROOT}/work/run_ourfm; cd $W
export TMPDIR=${PROJECT_ROOT}/work/tmp
L=$W/logs/finish.log
say() { echo "[finish] $(date '+%m-%d %H:%M') $*" | tee -a $L; }

say "START"

D=$W/fm/tchl/tchl_chr11_r2
if [ -s $D.ld.vars ] && [ ! -s $D.summary.tsv ]; then
  NLD=$(wc -l < $D.ld.vars); NST=$(( $(wc -l < $D.step2.txt) - 1 ))
  if [ "$NLD" -gt "$NST" ]; then
    mkdir -p $W/fm/tchl/stale
    for x in ld.bin ld.vars ld.json ld.done ld.py.log; do [ -e $D.$x ] && mv $D.$x $W/fm/tchl/stale/tchl_chr11_r2.$x; done
    say "chr11_r2 구버전 LD 격리 (ld=$NLD > step2=$NST) -> fm/tchl/stale/"
  fi
fi

run_region() {
  local RID=$1 T CH ROW S E
  T=${RID%%_*}; CH=$(echo $RID | sed 's/.*_chr\([0-9]*\)_.*/\1/')
  ROW=$(awk -F'\t' -v r="$RID" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next} $h["region_id"]==r{print $h["start"]"\t"$h["end"]; exit}' $W/regions/$T.chr$CH.regions.tsv)
  S=$(echo "$ROW" | cut -f1); E=$(echo "$ROW" | cut -f2)
  [ -n "$S" ] || { say "$RID 좌표 못 찾음 — 건너뜀"; return 1; }
  say "$RID 시작 (chr$CH $S-$E)"
  bash $W/scripts/fm_region.sh $T $CH $S $E $RID >> $W/logs/finish.$RID.log 2>&1
  local RC=$?
  say "$RID rc=$RC $(tail -1 $W/logs/finish.$RID.log | cut -c1-100)"
  return 0
}

for RID in dm_chr3_r1 tchl_chr11_r2 tchl_chr11_r3; do
  [ -s $W/fm/${RID%%_*}/$RID.summary.tsv ] && { say "$RID 이미 완료"; continue; }
  run_region $RID
done
say "STAGE1_DONE 완료구역=$(ls $W/fm/*/*.summary.tsv 2>/dev/null | grep -vc mi1000)"

R=/usr/bin/Rscript
for RID in dm_chr10_r2 lip_chr11_r1 tchl_chr19_r1 tchl_chr19_r2 tchl_chr20_r1; do
  T=${RID%%_*}; CH=$(echo $RID | sed 's/.*_chr\([0-9]*\)_.*/\1/'); D=$W/fm/$T/$RID
  [ -s $D.mi1000.summary.tsv ] && { say "$RID 재적합 이미 있음"; continue; }
  [ -s $D.step2.txt ] && [ -s $D.ld.bin ] || { say "$RID 재적합 불가(입력 없음)"; continue; }
  say "$RID 재적합 시작 (max_iter=1000)"
  FM_MAX_ITER=1000 $R $W/scripts/fm_susie.R $T $CH $RID $D.step2.txt $D.ld $D.mi1000 \
    >> $W/logs/finish.$RID.mi1000.log 2>&1
  say "$RID 재적합 rc=$? :: $(tail -1 $W/logs/finish.$RID.mi1000.log | cut -c1-110)"
done
say "STAGE2_DONE mi1000=$(ls $W/fm/*/*.mi1000.summary.tsv 2>/dev/null | wc -l)"

RID=htn_chr12_r2; T=htn; CH=12; D=$W/fm/$T/$RID
if [ ! -s $D.summary.tsv ]; then
  rm -f $D.skip
  ROW=$(awk -F'\t' -v r="$RID" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next} $h["region_id"]==r{print $h["start"]"\t"$h["end"]; exit}' $W/regions/$T.chr$CH.regions.tsv)
  S=$(echo "$ROW" | cut -f1); E=$(echo "$ROW" | cut -f2)
  say "$RID 단독 시작 (상한 해제, 예산 120G, chr$CH $S-$E)"
  MAXVAR=200000 FM_MAXVAR_MEM=200000 bash $W/scripts/fm_region.sh $T $CH $S $E $RID \
    >> $W/logs/finish.$RID.log 2>&1
  say "$RID rc=$? :: $(tail -1 $W/logs/finish.$RID.log | cut -c1-110)"
fi

python3 - <<'EOF' >> $L 2>&1
import os as _env_os, re as _env_re
def _env_path(value):
    return _env_re.sub(r"\$\{([A-Z][A-Z0-9_]*)\}", lambda match: _env_os.environ[match.group(1)], value)
import glob, csv
rows=[]
for f in sorted(glob.glob(_env_path('${PROJECT_ROOT}/work/run_ourfm/fm/*/*.summary.tsv'))):
    if 'mi1000' in f: continue
    r=list(csv.DictReader(open(f), delimiter='\t'))[0]; rows.append(r)
print('  최종 구역', len(rows), '| CS 합계', sum(int(r['n_cs']) for r in rows),
      '| 미수렴', sum(1 for r in rows if r['converged']!='TRUE'))
EOF
say "FINISH_DONE"
echo ok > $W/logs/finish.done
