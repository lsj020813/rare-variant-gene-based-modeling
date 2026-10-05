#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
T=${1:-htn}; N=22; NG=${2:-20}
W=${PROJECT_ROOT}/work; O=$W/ref/wtest; mkdir -p $O
S1=$W/ref/saige_step1_v4; GF=$W/ref/groupfiles_bwg/chr$N.B_3kb_re2g.txt
V=$W/ref/band_vcf/chr$N.band.vcf.gz; TRUTH=$W/ref/saige_step2_bwg/$T.chr$N
BCF=${BCFTOOLS:-bcftools}
L=$O/wtest.log; say(){ echo "[wtest] $(date '+%m-%d %H:%M') $*" | tee -a $L; }

for f in $S1/${T}_v4.rda $S1/${T}_v4.varianceRatio.txt $GF $V $TRUTH; do
  [ -s "$f" ] || { say "GATE FAIL: 없음 $f"; exit 2; }
done
UD=$(grep -m1 '^UD=' $W/run_bwg2/assoc.sh | cut -d= -f2- | tr -d '"')
[ -n "$UD" ] || { say "GATE FAIL: udocker 정의 못 읽음"; exit 3; }
say "START trait=$T genes=$NG UD=$(echo $UD | cut -c1-60)"

if [ ! -s $O/maf.tsv ]; then
  $BCF query -f '%CHROM:%POS:%REF:%ALT\t%INFO/MAF\n' $V > $O/maf.tsv 2>$O/maf.err
  say "MAF 표 $(wc -l < $O/maf.tsv) 변이"
fi
[ -s $O/maf.tsv ] || { say "GATE FAIL: MAF 표 비어있음 $(tail -1 $O/maf.err)"; exit 4; }

python3 - "$GF" "$O" "$NG" <<'PYEOF'
import sys
GF, O, NG = (sys.argv[1], sys.argv[2], int(sys.argv[3]))
maf = {}
for line in open(f'{O}/maf.tsv'):
    k, v = line.rstrip('\n').split('\t')
    try:
        maf[k] = float(v)
    except ValueError:
        pass
genes, order = ({}, [])
for line in open(GF):
    t = line.split()
    if len(t) < 3:
        continue
    g, kind = (t[0], t[1])
    if g not in genes:
        genes[g] = {}
        order.append(g)
    genes[g][kind] = t[2:]
sel = order[:NG]
beta = lambda m: 25.0 * (1.0 - m) ** 24  # dbeta(m; 1, 25)
n_missing = 0
with open(f'{O}/gf_A0.txt', 'w') as f0, open(f'{O}/gf_A1.txt', 'w') as f1, open(f'{O}/gf_A2.txt', 'w') as f2:
    for g in sel:
        v, a = (genes[g]['var'], genes[g].get('anno', ['all'] * len(genes[g]['var'])))
        assert len(v) == len(a), (g, len(v), len(a))
        for fh in (f0, f1, f2):
            fh.write(f'{g} var ' + ' '.join(v) + '\n')
            fh.write(f'{g} anno ' + ' '.join(a) + '\n')
        f1.write(f'{g} weight ' + ' '.join(['1.0'] * len(v)) + '\n')
        ws = []
        for k in v:
            m = maf.get(k)
            if m is None:
                n_missing += 1
                m = 0.005
            ws.append(f'{beta(m):.6g}')
        f2.write(f'{g} weight ' + ' '.join(ws) + '\n')
print(f"선택 유전자 {len(sel)} | 변이 합 {sum((len(genes[g]['var']) for g in sel))} | MAF 미검색 {n_missing}")
PYEOF
say "그룹파일 3종 작성 ($(awk '{print $1}' $O/gf_A0.txt | sort -u | wc -l) 유전자)"

run(){ local TAG=$1 G=$2 EXTRA=${3:-}
  local OUT=$O/$TAG
  [ -s $OUT ] && { say "$TAG 이미 있음"; return 0; }
  eval $UD step2_SPAtests.R \
    --vcfFile=\"$V\" --vcfFileIndex=\"$V.csi\" --vcfField=\"DS\" --chrom=\"$N\" --AlleleOrder=ref-first \
    --minMAF=0 --minMAC=0.5 --LOCO=FALSE --is_fastTest=TRUE \
    --GMMATmodelFile=\"$S1/${T}_v4.rda\" --varianceRatioFile=\"$S1/${T}_v4.varianceRatio.txt\" \
    --groupFile=\"$G\" --annotation_in_groupTest=\"all\" --maxMAF_in_groupTest=0.01 \
    --SAIGEOutputFile=\"$OUT\" $EXTRA > $OUT.log 2>&1
  local RC=$?; local NL=$(wc -l < $OUT 2>/dev/null || echo 0)
  say "$TAG rc=$RC rows=$NL $( [ "$NL" -lt 2 ] && tail -2 $OUT.log | tr '\n' ' ' | cut -c1-120 )"
}
run A0 $O/gf_A0.txt
run A1 $O/gf_A1.txt
run A2 $O/gf_A2.txt
run A3 $O/gf_A0.txt "--is_no_weight_in_groupTest=TRUE"

python3 - "$O" "$TRUTH" <<'PYEOF' | tee -a $L
import sys, csv, os
O, TRUTH = (sys.argv[1], sys.argv[2])

def load(f):
    d = {}
    if not (os.path.exists(f) and os.path.getsize(f) > 0):
        return d
    for r in csv.DictReader(open(f), delimiter='\t'):
        if r.get('Group') == 'all':
            try:
                d[r['Region']] = float(r['Pvalue'])
            except (ValueError, KeyError):
                pass
    return d
arms = {t: load(f'{O}/{t}') for t in ('A0', 'A1', 'A2', 'A3')}
truth = load(TRUTH)
print('  팔별 유전자 수:', {k: len(v) for k, v in arms.items()}, '| 정답지', len(truth))
keys = sorted(set(arms['A0']) & set(truth)) if arms['A0'] else []

def same(a, b, ks):
    return sum((1 for k in ks if abs(a[k] - b[k]) <= 1e-08 * max(1.0, abs(a[k]))))
if keys:
    print(f"  A0 vs 정답지: {same(arms['A0'], truth, keys)}/{len(keys)} 동일  ← 하네스 검증")
    for t in ('A1', 'A2', 'A3'):
        ks = sorted(set(arms['A0']) & set(arms[t]))
        if ks:
            print(f"  A0 vs {t}: {same(arms['A0'], arms[t], ks)}/{len(ks)} 동일")
    ks = sorted(set(arms.get('A1', {})) & set(arms.get('A3', {})))
    if ks:
        print(f"  A1 vs A3: {same(arms['A1'], arms['A3'], ks)}/{len(ks)} 동일")
    print(f"\n  {'유전자':24s}{'정답지':>12s}{'A0':>12s}{'A1':>12s}{'A2':>12s}{'A3':>12s}")
    for k in keys[:8]:
        row = f"  {k:24s}{truth[k]:12.4g}{arms['A0'][k]:12.4g}"
        for t in ('A1', 'A2', 'A3'):
            row += f"{arms[t].get(k, float('nan')):12.4g}"
        print(row)
    for t in ('A1', 'A2'):
        vals = [v for v in arms[t].values()]
        if vals:
            print(f'  {t} 병리 지표: P<1e-20 인 유전자 {sum((1 for v in vals if v < 1e-20))}/{len(vals)} · 최소 P {min(vals):.3g}')
PYEOF
say "WTEST_DONE"
echo ok > $O/wtest.done
