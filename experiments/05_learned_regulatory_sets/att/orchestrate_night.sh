#!/bin/bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -u
ATT=${PROJECT_ROOT}/work/fset/att
PY=python3
LOG=$ATT/logs/orchestrate.log
export TMPDIR=${PROJECT_ROOT}/work/tmp
mkdir -p $ATT/logs
say() { echo "$(date -Is) $*" | tee -a $LOG; }

say "ORCHESTRATE START"

for i in $(seq 1 360); do
  n=0
  for t in lip dm htn; do [ -f $ATT/trait_$t/att_fit.done ] && n=$((n+1)); done
  [ "$n" -ge 3 ] && break
  sleep 30
done
say "STEP1 형질 완료 $n/3"

for i in $(seq 1 720); do
  [ -f $ATT/b200/burden_stage.done ] && break
  sleep 30
done
if [ ! -f $ATT/b200/burden_stage.done ]; then
  say "STEP2 FAIL burden 미완료 -> B=200 fit 건너뜀"
else
  say "STEP2 burden 완료 ($(cat $ATT/b200/logs/burden_exit.log 2>/dev/null | wc -l)/22)"
  ncol=$(ls $ATT/b200/private/burden/col*.npy 2>/dev/null | wc -l)
  say "STEP2 consolidate 열 $ncol 개 (기대 212)"

  if [ "$ncol" -ge 212 ]; then
    export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
    cd $ATT
    for rng in "1 50" "51 100" "101 150" "151 200"; do
      set -- $rng
      say "STEP3 워커 시작 L=$1..$2"
      nice -n 19 ionice -c3 $PY att_fit_b200.py 200 $1 $2 \
        > $ATT/b200/logs/fit_$1_$2.log 2>&1 &
      sleep 20
    done
    wait
    say "STEP3 4워커 종료"

    $PY - <<'PYEOF' >> $LOG 2>&1
import os as _env_os, re as _env_re
def _env_path(value):
    return _env_re.sub(r"\$\{([A-Z][A-Z0-9_]*)\}", lambda match: _env_os.environ[match.group(1)], value)
import glob, csv, os
d = _env_path('${PROJECT_ROOT}/work/fset/att/b200')
parts = sorted(glob.glob(d + '/null_part_*.csv'))
rows = []
for p in parts:
    rows += list(csv.DictReader(open(p)))
if rows:
    rows.sort(key=lambda r: int(r['shuffle']))
    with open(d + '/att_label_shuffle_null.csv', 'w') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print('MERGE ok parts=%d rows=%d' % (len(parts), len(rows)))
else:
    print('MERGE FAIL: no parts')
PYEOF
    say "STEP4 병합 완료"
  else
    say "STEP3 SKIP: consolidate 열 부족($ncol)"
  fi
fi

cd $ATT
for d in b200 trait_lip trait_dm trait_htn lamfix; do
  if [ -f $ATT/$d/att_gene_results.csv ]; then
    $PY subset_bbj.py $d > $ATT/$d/subset_bbj.log 2>&1 \
      && say "STEP5 subset_bbj $d ok" || say "STEP5 subset_bbj $d FAIL"
  else
    say "STEP5 subset_bbj $d SKIP (결과 없음)"
  fi
done

$PY $ATT/night_report.py > $ATT/logs/report.log 2>&1 \
  && say "STEP6 보고서 생성 ok -> $ATT/MORNING_REPORT.md" \
  || say "STEP6 보고서 FAIL"

touch $ATT/night.stop
say "ORCHESTRATE DONE"
