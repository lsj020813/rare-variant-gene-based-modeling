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
import glob, gzip, collections, json
R = _configured('${PROJECT_ROOT}/work/ref')
THR = 2.5e-06
coord = {}
with gzip.open(f'{R}/deductive/gencode.sorted.gtf.gz', 'rt') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.split('\t')
        if len(f) < 9 or f[2] != 'gene' or f[0] not in ('19', 'chr19'):
            continue
        a = f[8]
        if 'gene_type "protein_coding"' not in a:
            continue
        gid = a.split('gene_id "')[1].split('"')[0].split('.')[0]
        coord[gid] = (int(f[3]), a.split('gene_name "')[1].split('"')[0] if 'gene_name "' in a else '?')

def sigset(files):
    out = {}
    for fp in files:
        with open(fp) as fh:
            h = fh.readline().rstrip('\n').split('\t')
            if 'Pvalue' not in h:
                continue
            gi = h.index('Region')
            pi = h.index('Pvalue')
            for line in fh:
                fl = line.rstrip('\n').split('\t')
                try:
                    p = float(fl[pi])
                except (ValueError, IndexError):
                    continue
                if p < THR:
                    g = fl[gi]
                    if g not in out or p < out[g]:
                        out[g] = p
    return out
S = {'A': sigset(sorted(glob.glob(f'{R}/saige_step2_v4/tchl.chr19.part[0-9][0-9][0-9]'))), 'B': sigset(sorted(glob.glob(f'{R}/saige_step2_pilot/tchl.chr19.B_3kb_re2g.B.part[0-9][0-9][0-9]'))), 'C': sigset([f'{R}/saige_step2_pilot/tchl.chr19.C_nearest']), 'D': sigset([f'{R}/saige_step2_pilot/tchl.chr19.D_3kb_only'])}

def loci(sig, win=500000):
    pts = sorted(((coord[g.split('.')[0]][0], g) for g in sig if g.split('.')[0] in coord))
    out = []
    cur = []
    for pos, g in pts:
        if cur and pos - cur[-1][0] > win:
            out.append(cur)
            cur = []
        cur.append((pos, g))
    if cur:
        out.append(cur)
    return out

def label(L):
    return '/'.join(sorted({coord[g.split('.')[0]][1] for _, g in L})[:3])
print('=== 팔별 유의 좌위 (거리 500kb) ===')
LOC = {}
for k, sig in S.items():
    ls = loci(sig)
    LOC[k] = [label(L) for L in ls]
    print(f'  {k}: {len(sig)} genes -> {len(ls)} loci : {LOC[k]}')
print()
print('=== 좌위 중심위치로 대조 (500kb 창) ===')

def centers(sig):
    return sorted((sum((p for p, _ in L)) / len(L) for L in loci(sig)))
for k in S:
    print(f"  {k}: {[f'{c / 1000000.0:.2f}Mb' for c in centers(S[k])]}")
print()
print('=== 고유 발견: B 에만 있고 C 에 없는 유전자 ===')
b_only = set(S['B']) - set(S['C'])
c_only = set(S['C']) - set(S['B'])
for g in sorted(b_only, key=lambda x: S['B'][x]):
    print(f"  B-only {coord.get(g.split('.')[0], (0, '?'))[1]:<12} p={S['B'][g]:.2e}")
for g in sorted(c_only, key=lambda x: S['C'][x]):
    print(f"  C-only {coord.get(g.split('.')[0], (0, '?'))[1]:<12} p={S['C'][g]:.2e}")
print()
print('=== B vs D (rE2G 순기여) ===')
bd = set(S['B']) - set(S['D'])
db = set(S['D']) - set(S['B'])
print(f"  B-only vs D: {[coord.get(g.split('.')[0], (0, '?'))[1] for g in bd]}")
print(f"  D-only vs B: {[coord.get(g.split('.')[0], (0, '?'))[1] for g in db]}")
json.dump({k: {g: p for g, p in v.items()} for k, v in S.items()}, open(_configured('${PROJECT_ROOT}/work/run_audit/arm_sig.json'), 'w'), indent=1)
print('\nwritten arm_sig.json')
