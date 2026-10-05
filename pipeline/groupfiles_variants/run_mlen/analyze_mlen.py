#!/usr/bin/env python3
import json, math, sys
import numpy as np
rows = json.load(open(sys.argv[1]))
m = np.array([r['m'] for r in rows], float)
L = np.array([r['L'] for r in rows], float)
mp = np.array([r['m_prom'] for r in rows], float)
sym = [r['sym'] for r in rows]
print(f'genes {len(rows):,}')
print('\n=== (1) m 분포 ===')
print(f'  median {np.median(m):.0f} | p90 {np.percentile(m, 90):.0f} | max {m.max():.0f} ({sym[int(np.argmax(m))]}) | mean {m.mean():.2f} | min {m.min():.0f}')
print(f'  max/median {m.max() / np.median(m):.1f}x | mean>median? {m.mean() > np.median(m)}')
print(f'  L: median {np.median(L) / 1000:.1f}kb | p90 {np.percentile(L, 90) / 1000:.1f}kb | max {L.max() / 1000000.0:.2f}Mb ({sym[int(np.argmax(L))]})')

def rankdata(a):
    o = np.argsort(a, kind='mergesort')
    r = np.empty(len(a), float)
    r[o] = np.arange(1, len(a) + 1)
    _, inv, cnt = np.unique(a, return_inverse=True, return_counts=True)
    sums = np.zeros(len(cnt))
    np.add.at(sums, inv, r)
    return (sums / cnt)[inv]
rm, rl = (rankdata(m), rankdata(L))
rho = np.corrcoef(rm, rl)[0, 1]
print(f'\n=== (2) Spearman rho(m, L) = {rho:.4f} ===')
dens = m / (L / 1000.0)
print('\n=== (3) 밀도 m/L (변이/kb) ===')
print(f'  median {np.median(dens):.3f} | IQR {np.percentile(dens, 25):.3f}~{np.percentile(dens, 75):.3f} | p5 {np.percentile(dens, 5):.3f} | p95 {np.percentile(dens, 95):.3f}')
print(f'  CV(변동계수) {dens.std() / dens.mean():.3f}  <- 0 에 가까우면 순수 길이 문제')
print(f'  p95/p5 {np.percentile(dens, 95) / max(np.percentile(dens, 5), 1e-09):.1f}x')
ok = (m > 0) & (L > 0)
lx, ly = (np.log10(L[ok]), np.log10(m[ok]))
slope, intercept = np.polyfit(lx, ly, 1)
pred = slope * lx + intercept
r2 = 1 - ((ly - pred) ** 2).sum() / ((ly - ly.mean()) ** 2).sum()
print(f'\n=== (4) log10 m ~ log10 L ===')
print(f'  기울기 {slope:.4f} | 절편 {intercept:.4f} | R2 {r2:.4f}  (기울기 1 = 비례)')
resid = ly - pred
idx = np.argsort(-resid)[:20]
syms_ok = [s for s, k in zip(sym, ok) if k]
Lok, mok = (L[ok], m[ok])
print('\n=== (5) 잔차 상위 20 (길이 대비 m 과다) ===')
print(f"  {'gene':<14}{'m':>7}{'L(kb)':>10}{'m/kb':>8}{'resid':>8}")
for i in idx:
    print(f'  {syms_ok[i]:<14}{mok[i]:>7.0f}{Lok[i] / 1000:>10.1f}{mok[i] / (Lok[i] / 1000):>8.2f}{resid[i]:>8.2f}')
frac = np.divide(mp, m, out=np.zeros_like(mp), where=m > 0)
print('\n=== (6) m 의 프로모터(TSS+-3kb) 비중 ===')
print(f'  전체 합계: 프로모터 {mp.sum():,.0f} / 전체 {m.sum():,.0f} = {mp.sum() / m.sum() * 100:.1f}%')
print(f'  유전자별 비중: median {np.median(frac) * 100:.1f}% | mean {frac.mean() * 100:.1f}% | 프로모터만인 유전자 {(frac >= 0.999).sum():,} | 프로모터 0인 유전자 {(frac <= 0.001).sum():,}')
json.dump({'genes': len(rows), 'm_median': float(np.median(m)), 'm_p90': float(np.percentile(m, 90)), 'm_max': float(m.max()), 'm_max_gene': sym[int(np.argmax(m))], 'm_mean': float(m.mean()), 'rho_m_L': float(rho), 'dens_median': float(np.median(dens)), 'dens_cv': float(dens.std() / dens.mean()), 'loglog_slope': float(slope), 'loglog_r2': float(r2), 'prom_frac_overall': float(mp.sum() / m.sum()), 'prom_frac_median': float(np.median(frac))}, open(sys.argv[1].replace('.json', '_summary.json'), 'w'), indent=1)
print('\nANALYZE_DONE')
