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
import gzip, glob, json, math, os, re, sys
from collections import defaultdict
GTF = _configured('${PROJECT_ROOT}/work/ref/deductive/gencode.sorted.gtf.gz')
GP = _configured('${PROJECT_ROOT}/work/ref/groupfiles_pilot')
BAND = _configured('${PROJECT_ROOT}/work/ref/band_vcf')
OUT = _configured('${PROJECT_ROOT}/work/run_mlen')
HALF = 3000
os.makedirs(OUT, exist_ok=True)
genes = {}
with gzip.open(GTF, 'rt') as fh:
    for line in fh:
        if line[0] == '#':
            continue
        f = line.split('\t', 9)
        if len(f) < 9 or f[2] != 'gene':
            continue
        if 'gene_type "protein_coding"' not in f[8]:
            continue
        i = f[8].find('gene_id "')
        gid = f[8][i + 9:f[8].find('"', i + 9)]
        j = f[8].find('gene_name "')
        sym = f[8][j + 11:f[8].find('"', j + 11)] if j >= 0 else gid
        s, e = (int(f[3]), int(f[4]))
        tss = s if f[6] == '+' else e
        genes[gid] = (f[0], tss, s, e, sym)
print(f'protein-coding genes: {len(genes):,}', flush=True)
assert genes, 'GATE FAIL: no genes'

def load_group(path):
    per = {}
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 3 or f[1] != 'var':
                continue
            per[f[0]] = set(f[2:])
        return per
files = sorted(glob.glob(f'{GP}/chr*.B_3kb_re2g.txt'))
assert files, f'GATE FAIL: no group files under {GP}'
print(f'group files: {[os.path.basename(x) for x in files]}', flush=True)
rows = []
for path in files:
    per = load_group(path)
    for gid, keys in per.items():
        gi = genes.get(gid) or genes.get(next((g for g in genes if g.split('.')[0] == gid.split('.')[0]), ''), None)
        if not gi:
            continue
        _, tss, gs, ge, sym = gi
        lo = min(tss - HALF, gs)
        hi = max(tss + HALF, ge)
        L = hi - lo + 1
        pl, ph = (tss - HALF, tss + HALF)
        m_prom = 0
        for k in keys:
            p = int(k.split(':')[1])
            if pl <= p <= ph:
                m_prom += 1
        rows.append((gid, sym, len(keys), L, m_prom, len(keys) - m_prom))
assert rows, 'GATE FAIL: no gene rows matched between group files and GTF'
json.dump([{'gid': g, 'sym': s, 'm': m, 'L': L, 'm_prom': mp, 'm_body': mb} for g, s, m, L, mp, mb in rows], open(f'{OUT}/mlen_rows.json', 'w'))
print(f'MLEN_ROWS {len(rows)}', flush=True)
print('MLEN_PROBE_DONE', flush=True)
