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
import glob, gzip, collections
R = _configured('${PROJECT_ROOT}/work/ref')
GID = 'ENSG00000129354'
TSS = _required_int('TSS_POSITION')

def gv(fp):
    for line in open(fp):
        f = line.split()
        if len(f) > 2 and f[1] == 'var' and (f[0].split('.')[0] == GID):
            return set(f[2:])
    return set()
B = gv(f'{R}/groupfiles_pilot/chr19.B_3kb_re2g.txt')
D = gv(f'{R}/groupfiles_pilot/chr19.D_3kb_only.txt')
extra = sorted(B - D, key=lambda k: int(k.split(':')[1]))
print(f'B={len(B)}  D={len(D)}  rE2G-only={len(extra)}')
print('\n=== rE2G-only variant positions (relative to configured TSS) ===')
for k in extra:
    pos = int(k.split(':')[1])
    print(f'  {k:<28} pos {pos:,}  d={pos - TSS:+,}')
sa = sorted(glob.glob(f'{R}/saige_step2_pilot/tchl.chr19.B_3kb_re2g.B.part*.singleAssoc.txt'))
print(f'\nsingleAssoc files: {len(sa)}')
rows = []
for fp in sa:
    with open(fp) as fh:
        h = fh.readline().rstrip('\n').split('\t')
        if 'p.value' not in h and 'Pvalue' not in h:
            print('  cols:', h[:12])
            continue
        pi = h.index('p.value') if 'p.value' in h else h.index('Pvalue')
        mi = h.index('MarkerID') if 'MarkerID' in h else 0
        af = h.index('AF_Allele2') if 'AF_Allele2' in h else None
        for line in fh:
            fl = line.rstrip('\n').split('\t')
            if fl[mi] in B:
                try:
                    rows.append((float(fl[pi]), fl[mi], fl[af] if af is not None else ''))
                except (ValueError, IndexError):
                    pass
rows.sort()
print(f"\n=== AP1M2's own variants, single-variant p (top 12 of {len(rows)}) ===")
for p, m, a in rows[:12]:
    tag = 'rE2G-only' if m in B - D else 'distance'
    print(f'  {m:<28} p={p:.2e}  AF={a:<10} [{tag}]')
n_rare_sig = sum((1 for p, _, _ in rows if p < 0.05))
print(f'\n  variants with single-variant p<0.05: {n_rare_sig}/{len(rows)}')
print(f'  of the rE2G-only ones: {sum((1 for p, m, _ in rows if p < 0.05 and m in B - D))}/{len(B - D)}')
