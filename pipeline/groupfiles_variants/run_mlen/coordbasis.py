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
import gzip, os, subprocess
import numpy as np
R = _configured('${PROJECT_ROOT}/work/ref')
BCF = _configured('${BCFTOOLS}')
HALF = 3000
genes = {}
with gzip.open(f'{R}/deductive/gencode.sorted.gtf.gz', 'rt') as fh:
    for line in fh:
        if line[0] == '#':
            continue
        f = line.split('\t', 9)
        if len(f) < 9 or f[2] != 'gene' or 'gene_type "protein_coding"' not in f[8]:
            continue
        i = f[8].find('gene_id "')
        gid = f[8][i + 9:f[8].find('"', i + 9)]
        s, e = (int(f[3]), int(f[4]))
        genes[gid] = (s if f[6] == '+' else e, s, e)
for N in (20, 13, 19):
    per = {}
    want = set()
    for line in open(f'{R}/groupfiles_bwg/chr{N}.B_3kb_re2g.txt'):
        f = line.split()
        if len(f) < 3 or f[1] != 'var':
            continue
        per[f[0]] = f[2:]
        want.update(f[2:])
    pos38 = {}
    p = subprocess.Popen(f"{BCF} query -f '%ID\t%POS\n' {R}/lift38_keyed/chr{N}.keyed38.vcf.gz", shell=True, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for line in p.stdout:
        k, v = line.rstrip('\n').split('\t', 1)
        kk = k[3:] if k.startswith('chr') else k
        if kk in want:
            pos38[kk] = int(v)
    p.wait()
    in38 = out38 = in37 = 0
    for gid, keys in per.items():
        gi = genes.get(gid)
        if not gi:
            continue
        tss, gs, ge = gi
        lo, hi = (min(tss - HALF, gs), max(tss + HALF, ge))
        for k in keys:
            p38 = pos38.get(k)
            p37 = int(k.split(':')[1])
            if p38 is not None:
                if lo <= p38 <= hi:
                    in38 += 1
                else:
                    out38 += 1
            if lo <= p37 <= hi:
                in37 += 1
    tot = in38 + out38
    print(f'chr{N}: 매핑 {tot:,} | 38좌표가 창 안 {in38:,} ({in38 / tot * 100:.1f}%) | 창 밖 {out38:,} ({out38 / tot * 100:.1f}%) || 37좌표를 창에 대면 {in37 / tot * 100:.1f}%')
print('COORDBASIS_DONE')
