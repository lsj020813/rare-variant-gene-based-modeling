import os as _config_os
import re as _config_re

def _required(name):
    value = _config_os.environ.get(name)
    if value is None or not value.strip():
        raise RuntimeError(f'Set {name} before running this script')
    return value

def _required_int(name):
    value = int(_required(name))
    if value <= 0:
        raise ValueError(f'{name} must be a positive integer')
    return value

def _configured(value):
    defaults = {'BCFTOOLS': 'bcftools', 'PLINK2': 'plink2', 'PYTHON': 'python3', 'TABIX': 'tabix', 'BGZIP': 'bgzip', 'SAMTOOLS': 'samtools', 'BEDTOOLS': 'bedtools', 'CROSSMAP': 'CrossMap', 'UDOCKER_BIN': 'udocker'}
    def resolve(match):
        name = match.group(1) or match.group(2)
        if name in _config_os.environ:
            return _required(name)
        return defaults.get(name) or _required(name)
    return _config_re.sub(r"\$\{([A-Z][A-Z0-9_]*)\}|\$([A-Z][A-Z0-9_]*)", resolve, value)

def _source_path(base, value):
    return value if _config_os.path.isabs(value) else _config_os.path.join(base, value)
import gzip, json, os, subprocess
import numpy as np
GTF = _configured('${PROJECT_ROOT}/work/ref/deductive/gencode.sorted.gtf.gz')
OUT = _configured('${PROJECT_ROOT}/work/run_mlen')
hdr = []
with gzip.open(GTF, 'rt') as fh:
    for i, l in enumerate(fh):
        if l[0] != '#':
            break
        hdr.append(l.strip())
print('=== (1) GTF 헤더 ===')
for h in hdr[:8]:
    print('  ', h[:110])
genes = {}
with gzip.open(GTF, 'rt') as fh:
    for line in fh:
        if line[0] == '#':
            continue
        f = line.split('\t', 9)
        if len(f) < 9 or f[2] != 'gene' or 'gene_type "protein_coding"' not in f[8]:
            continue
        i = f[8].find('gene_id "')
        gid = f[8][i + 9:f[8].find('"', i + 9)]
        j = f[8].find('gene_name "')
        sym = f[8][j + 11:f[8].find('"', j + 11)] if j >= 0 else gid
        s, e = (int(f[3]), int(f[4]))
        tss = s if f[6] == '+' else e
        genes[gid] = (f[0], tss, s, e, sym)
for gid, (c_, t, s, e, sym) in genes.items():
    if sym == 'TP53':
        print(f'  TP53 {c_}:{s:,}-{e:,}  -> ' + ('GRCh38' if s > 7600000 else 'GRCh37'))
        break

def Lstats(half):
    L = np.array([max(t + half, e) - min(t - half, s) + 1 for _, t, s, e, _ in genes.values()], float)
    return L
for half, lab in ((100000, '정답지 A (±100kb∪본체)'), (3000, '신 규칙 B (±3kb∪본체)')):
    L = Lstats(half)
    print(f'\n=== (2) {lab} ===')
    print(f'  L 중앙 {np.median(L) / 1000:.1f}kb | p10 {np.percentile(L, 10) / 1000:.1f}kb | p90 {np.percentile(L, 90) / 1000:.1f}kb | max {L.max() / 1000000.0:.2f}Mb | 합계 {L.sum() / 1000000000.0:.3f}Gb')
rows = json.load(open(f'{OUT}/mlen_rows_wg.json'))
m = np.array([r['m'] for r in rows], float)
L = np.array([r['L'] for r in rows], float)
sym = [r['sym'] for r in rows]
ok = (m > 0) & (L > 0)
lx, ly = (np.log10(L[ok]), np.log10(m[ok]))
sl, ic = np.polyfit(lx, ly, 1)
pred = 10 ** (sl * lx + ic)
ratio = m[ok] / pred
syms = [s for s, k in zip(sym, ok) if k]
Lo, mo = (L[ok], m[ok])
idx = np.argsort(-ratio)[:20]
print(f'\n=== (3) 잔차 재순위 (m / 예측), 기울기 {sl:.4f} ===')
print(f"  {'gene':<14}{'m':>7}{'예측':>8}{'배수':>7}{'L(kb)':>9}")
for i in idx:
    print(f'  {syms[i]:<14}{mo[i]:>7.0f}{pred[i]:>8.0f}{ratio[i]:>7.2f}{Lo[i] / 1000:>9.1f}')
csmd = [i for i, s in enumerate(syms) if s == 'CSMD1']
if csmd:
    i = csmd[0]
    r_ = int((ratio > ratio[i]).sum()) + 1
    print(f'  [참고] CSMD1: m {mo[i]:.0f} 예측 {pred[i]:.0f} 배수 {ratio[i]:.2f} 순위 {r_}/{len(ratio)}')
json.dump({'syms_top20': [syms[i] for i in idx], 'ratio_top20': [float(ratio[i]) for i in idx], 'slope': float(sl)}, open(f'{OUT}/resid_rerank.json', 'w'), indent=1)
print('PREPILOT_DONE')
