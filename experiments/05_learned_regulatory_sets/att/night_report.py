import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import csv
import glob
import json
import os
import numpy as np

ATT = _config_path('${PROJECT_ROOT}/work/fset/att')
OUT = ATT + '/MORNING_REPORT.md'
L = []

def w(s=''):
    L.append(s)

def gene_counts(p):
    if not os.path.exists(p):
        return None
    g = list(csv.DictReader(open(p)))
    arms = ['FLAT', 'FIXED', 'SIZE', 'ATT', 'INT', 'ATT_GENEHOLD',
            'ATT_MAF001']
    out = {a: sum(1 for r in g if _f(r.get('q_' + a)) < 0.05) for a in arms}
    out['_n'] = len(g)
    out['_d_att'] = float(np.nanmean([_f(r['delta_att_flat']) for r in g]))
    out['_d_int'] = float(np.nanmean([_f(r['delta_int_flat']) for r in g]))
    return out

def _f(x):
    try:
        v = float(x)
        return v if np.isfinite(v) else np.nan
    except (TypeError, ValueError):
        return np.nan

def null_summary(p):
    if not os.path.exists(p):
        return None
    r = list(csv.DictReader(open(p)))
    ns = np.array([int(x['n_sig_fdr05']) for x in r])
    ni = np.array([int(x['n_sig_int_fdr05']) for x in r])
    md = np.array([float(x['mean_delta']) for x in r])
    wm = np.array([float(x['mean_w_major']) for x in r])
    return dict(B=len(r), att_mean=ns.mean(), att_rng=(ns.min(), ns.max()),
                int_mean=ni.mean(), int_rng=(ni.min(), ni.max()),
                md_mean=md.mean(), w_mean=wm.mean(),
                floor=1.0 / (len(r) + 1))

def wstats(p):
    if not os.path.exists(p):
        return None
    v = np.array([_f(r['w_major_att'])
                  for r in csv.DictReader(open(p))])
    v = v[np.isfinite(v)]
    return dict(n=len(v), mean=v.mean(), sd=v.std(), lo=v.min(), hi=v.max(),
                f05=float((np.abs(v - .5) > .05).mean()),
                f20=float((np.abs(v - .5) > .2).mean()))

def pv(obs, null):
    null = np.asarray(null, float)
    return (1 + int((null >= obs).sum())) / (len(null) + 1)

w('# ATT 재검정 — 야간 자동 보고서')
w()
w('생성: `%s`' % __import__('time').strftime('%Y-%m-%d %H:%M:%S'))
w()
w('## 0. 한 줄 요약')
w()
w('원본 `fset/att` 의 ATT 결과는 λ 그리드 오설정으로 게이트가 죽어 있어(전 유전자 w=0.5) '
  '**미검정 상태였다.** λ 를 목적함수 스케일에 맞춰 고친 뒤 재검정했고, 양성 대조로 '
  '검정력이 있음을 확인했다. 아래가 그 결과다.')
w()

w('## 1. 게이트가 살아났는가 (w 분포)')
w()
w('| 실행 | λ 그리드 | w 평균 | w 표준편차 | 범위 | \\|w-.5\\|>.05 | \\|w-.5\\|>.2 |')
w('|---|---|---|---|---|---|---|')
for tag, d, lam in [('원본(1-SE)', 'att', '[1e-3..1.0]'),
                    ('CV-min 1차', 'cvmin', '[1e-5..1.0]'),
                    ('**λ수정**', 'lamfix', '**[0..1e-5]**'),
                    ('B=200', 'b200', '[0..1e-5]')]:
    s = wstats(ATT + '/att_weights.csv' if d == 'att'
               else '%s/%s/att_weights.csv' % (ATT, d))
    if s:
        w('| %s | %s | %.4f | %.4f | %.4f~%.4f | %.1f%% | %.1f%% |'
          % (tag, lam, s['mean'], s['sd'], s['lo'], s['hi'],
             100 * s['f05'], 100 * s['f20']))
w()
w('`probe_lambda_scale.py` 실측: 목적함수 evidence 는 ~2.5e-05 스케일인데 벌점은 '
  'λ·‖params‖² 이고 ‖params‖² 는 O(1)~O(100). **λ≥1e-5 면 파라미터가 0 으로 눌려 '
  'gate_w=sigmoid(0)=0.5 로 고정**된다. 원본 그리드는 전 구간이 이 죽은 영역이었다.')
w()

w('## 2. 양성 대조 — 검정력이 있는가')
w()
rows = []
for f in (ATT + '/poscontrol/pos_control.csv',
          ATT + '/poscontrol_hi/pos_control.csv'):
    if os.path.exists(f):
        rows += list(csv.DictReader(open(f)))
if rows:
    for scen in ('P1_MAJOR', 'P2_CCRE'):
        w('**%s** (%s)' % (scen, '신호=항상 major, 최적 규칙은 상수'
                           if scen == 'P1_MAJOR'
                           else '신호=frac_ccre 큰 쪽, 규칙 학습 필요'))
        w()
        w('| δ | w_causal | toward(정답방향) | ATT 유의 | FLAT 유의 | dz(인과) |')
        w('|---|---|---|---|---|---|')
        rr = sorted([r for r in rows if r['scenario'] == scen],
                    key=lambda r: (float(r['delta']), int(r['rep'])))
        for r in rr:
            w('| %s (rep%s) | %.3f | %.2f | %s | %s | %+.3f |'
              % (r['delta'], r['rep'], float(r['w_causal_mean']),
                 float(r['frac_w_toward_truth']), r['n_sig_att'],
                 r['n_sig_flat'], float(r['mean_delta_z_causal'])))
        w()
    w('**해석**: δ=0 에서 교정 정상(유의 0, toward 우연 수준). δ≥0.1 에서 게이트가 '
      '정답 규칙을 거의 완벽히 학습(toward→1.00)하고, 특징 의존 규칙(P2)에서 '
      'ATT 가 FLAT 을 크게 앞선다. **검출 하한 δ≈0.1.**')
else:
    w('_양성 대조 산출물 없음_')
w()

w('## 3. 본 검정 — 형질별 결과')
w()
w('| 실행 | 형질 | N | FLAT | ATT | INT | ATT 귀무평균[범위] | 순열 p(ATT n_sig) |')
w('|---|---|---|---|---|---|---|---|')
for d, tr in [('lamfix', 'TCHL(양적)'), ('b200', 'TCHL(양적) B=200'),
              ('trait_lip', 'LIP(이진)'), ('trait_dm', 'DM(이진)'),
              ('trait_htn', 'HTN(이진)')]:
    gc = gene_counts('%s/%s/att_gene_results.csv' % (ATT, d))
    ns = null_summary('%s/%s/att_label_shuffle_null.csv' % (ATT, d))
    if not gc:
        continue
    if ns:
        r = list(csv.DictReader(
            open('%s/%s/att_label_shuffle_null.csv' % (ATT, d))))
        p = pv(gc['ATT'], [int(x['n_sig_fdr05']) for x in r])
        w('| %s | %s | %d | %d | %d | %d | %.1f [%d, %d] | %.4f (B=%d) |'
          % (d, tr, gc['_n'], gc['FLAT'], gc['ATT'], gc['INT'],
             ns['att_mean'], ns['att_rng'][0], ns['att_rng'][1], p, ns['B']))
    else:
        w('| %s | %s | %d | %d | %d | %d | _귀무 없음_ | — |'
          % (d, tr, gc['_n'], gc['FLAT'], gc['ATT'], gc['INT']))
w()

w('## 4. BBJ 부분집합 (외부 증거로 좁힌 22개 도메인)')
w()
w('선택 근거가 외부 코호트(BBJ TC FINEMAP PIP≥0.1, ±500kb)라 우리 표현형 잡음과 '
  '무관 → 순환논증 없음. FDR 은 부분집합 안에서 재계산.')
w()
w('| 실행 | 범위 | 통계량 | 관측 | 귀무평균 | 귀무범위 | 순열 p |')
w('|---|---|---|---|---|---|---|')
for d in ('lamfix', 'b200', 'trait_lip', 'trait_dm', 'trait_htn'):
    p = '%s/%s/subset_bbj.json' % (ATT, d)
    if not os.path.exists(p):
        continue
    j = json.load(open(p))
    for scope, dd in j['result'].items():
        for k in ('ATT_n_sig', 'ATT_mean_delta'):
            v = dd['stats'][k]
            w('| %s | %s | %s | %.4f | %.4f | [%.3f, %.3f] | **%.4f** |'
              % (d, scope, k, v['obs'], v['null_mean'], v['null_min'],
                 v['null_max'], v['p']))
w()

w('## 5. 야간 자원 사용')
w()
mp = ATT + '/logs/night_monitor.log'
if os.path.exists(mp):
    ls = [x for x in open(mp) if 'load=' in x]
    if ls:
        def grab(line, key):
            for tok in line.split():
                if tok.startswith(key):
                    return tok.split('=', 1)[1]
            return ''
        av = [float(grab(x, 'availGB=')) for x in ls if grab(x, 'availGB=')]
        rs = [float(grab(x, 'ourRSS=').rstrip('G')) for x in ls
              if grab(x, 'ourRSS=')]
        l1 = [float(x.split('load=')[1].split()[0]) for x in ls]
        st = [int(grab(x, 'stopped=') or 0) for x in ls]
        w('- 표본 %d개 (5분 간격, %s ~ %s)'
          % (len(ls), ls[0][:19], ls[-1][:19]))
        w('- load1: 최소 %.2f / 중앙 %.2f / **최대 %.2f**'
          % (min(l1), float(np.median(l1)), max(l1)))
        w('- 가용 RAM: **최소 %.0fG** / 중앙 %.0fG' % (min(av), np.median(av)))
        w('- 우리 RSS 합: 중앙 %.1fG / **최대 %.1fG**' % (np.median(rs), max(rs)))
        w('- 워치독에 의한 일시정지(SIGSTOP) 발생 틱: **%d회**' % sum(1 for x in st if x))
w()
for wd, nm in [('logs/wd_fit.log', 'fit'), ('logs/wd_trait.log', 'trait'),
               ('b200/logs/watchdog.log', 'burden')]:
    p = ATT + '/' + wd
    if os.path.exists(p):
        k = sum(1 for x in open(p) if 'KILL' in x)
        s = sum(1 for x in open(p) if ' STOP ' in x)
        w('- 워치독 `%s`: STOP %d회, KILL %d회' % (nm, s, k))
w()
w('## 6. 산출물 위치')
w()
for d, desc in [('att', '원본(1-SE, 무효 — README 참조)'),
                ('cvmin', 'CV-min 1차(무효 확인용)'),
                ('lamfix', 'λ 수정 20셔플 ★주 결과'),
                ('b200', 'λ 수정 200셔플'),
                ('poscontrol', '양성 대조 저효과'),
                ('poscontrol_hi', '양성 대조 고효과'),
                ('trait_lip', 'LIP'), ('trait_dm', 'DM'),
                ('trait_htn', 'HTN')]:
    p = ATT if d == 'att' else ATT + '/' + d
    if os.path.exists(p):
        w('- `fset/att%s` — %s' % ('' if d == 'att' else '/' + d, desc))
w()

open(OUT, 'w').write('\n'.join(L) + '\n')
print('wrote', OUT, len(L), 'lines')
