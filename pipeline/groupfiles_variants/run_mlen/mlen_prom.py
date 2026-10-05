#!/usr/bin/env python3
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
import gzip, glob, json, os, subprocess, sys
import numpy as np
GTF = _configured('${PROJECT_ROOT}/work/ref/deductive/gencode.sorted.gtf.gz')
GP = _configured('${PROJECT_ROOT}/work/ref/groupfiles_bwg')
KEYED = _configured('${PROJECT_ROOT}/work/ref/lift38_keyed')
BCF = _configured('${BCFTOOLS}')
OUT = _configured('${PROJECT_ROOT}/work/run_mlen')
HALF = 3000
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
        genes[gid] = (f[0], s if f[6] == '+' else e, s, e, sym)
print(f'genes {len(genes):,}', flush=True)
rows = []
for N in range(1, 23):
    gp = f'{GP}/chr{N}.B_3kb_re2g.txt'
    kf = f'{KEYED}/chr{N}.keyed38.vcf.gz'
    if not (os.path.exists(gp) and os.path.exists(kf)):
        print(f'  chr{N}: SKIP (missing input)', flush=True)
        continue
    per, want = ({}, set())
    with open(gp) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 3 or f[1] != 'var':
                continue
            per[f[0]] = f[2:]
            want.update(f[2:])
    pos38 = {}
    p = subprocess.Popen(f"{BCF} query -f '%ID\\t%POS\\n' {kf}", shell=True, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for line in p.stdout:
        k, v = line.rstrip('\n').split('\t', 1)
        kk = k[3:] if k.startswith('chr') else k
        if kk in want:
            pos38[kk] = int(v)
    p.wait()
    hit = len(pos38)
    assert hit > 0, f'GATE FAIL chr{N}: no 37->38 key matched ({len(want):,} wanted)'
    n_gene = 0
    for gid, keys in per.items():
        gi = genes.get(gid)
        if not gi:
            continue
        _, tss, gs, ge, sym = gi
        lo, hi = (min(tss - HALF, gs), max(tss + HALF, ge))
        pl, ph = (tss - HALF, tss + HALF)
        m_prom = m_map = 0
        for k in keys:
            p38 = pos38.get(k)
            if p38 is None:
                continue
            m_map += 1
            if pl <= p38 <= ph:
                m_prom += 1
        rows.append({'gid': gid, 'sym': sym, 'm': len(keys), 'm_mapped': m_map, 'L': hi - lo + 1, 'm_prom': m_prom, 'm_body': m_map - m_prom, 'chrom': N})
        n_gene += 1
    print(f'  chr{N}: genes {n_gene:,} | keys wanted {len(want):,} | mapped {hit:,} ({hit / len(want) * 100:.1f}%)', flush=True)
assert rows, 'GATE FAIL: no rows'
json.dump(rows, open(f'{OUT}/mlen_prom_rows.json', 'w'))
m = np.array([r['m_mapped'] for r in rows], float)
mp = np.array([r['m_prom'] for r in rows], float)
frac = np.divide(mp, m, out=np.zeros_like(mp), where=m > 0)
print(f'\n=== (6) 정정판 — 프로모터(TSS+-3kb) 비중 ===')
print(f'  대조 가능 유전자 {len(rows):,} | 좌표 매핑된 변이 {m.sum():,.0f}')
print(f'  전체: 프로모터 {mp.sum():,.0f} / {m.sum():,.0f} = {mp.sum() / m.sum() * 100:.2f}%')
print(f'  유전자별 비중: median {np.median(frac) * 100:.2f}% | mean {frac.mean() * 100:.2f}% | p90 {np.percentile(frac, 90) * 100:.2f}%')
print(f'  프로모터 변이 0개인 유전자 {(mp == 0).sum():,} ({(mp == 0).mean() * 100:.1f}%)')
print(f'  프로모터만인 유전자 {(frac >= 0.999).sum():,}')
print('MLEN_PROM_DONE', flush=True)
