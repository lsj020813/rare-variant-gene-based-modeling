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
TARGET = 'ENSG00000129354'
gid = None
g_start = g_end = None
with gzip.open(f'{R}/deductive/gencode.sorted.gtf.gz', 'rt') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.split('\t')
        if len(f) < 9 or f[2] != 'gene' or f[0] not in ('19', 'chr19'):
            continue
        if 'gene_name "AP1M2"' in f[8]:
            gid = f[8].split('gene_id "')[1].split('"')[0]
            g_start, g_end, strand = (int(f[3]), int(f[4]), f[6])
            break
print(f'AP1M2 gene_id={gid}  {g_start:,}-{g_end:,} ({strand})')

def gv(fp, want):
    for line in open(fp):
        f = line.split()
        if len(f) > 2 and f[1] == 'var' and (f[0] == want):
            return set(f[2:])
    return set()
B = gv(f'{R}/groupfiles_pilot/chr19.B_3kb_re2g.txt', gid)
D = gv(f'{R}/groupfiles_pilot/chr19.D_3kb_only.txt', gid)
print(f'arm B variants: {len(B)}   arm D variants: {len(D)}   B-only: {len(B - D)}')
tss = g_start if strand == '+' else g_end
dists = sorted((abs(int(k.split(':')[1]) - tss) for k in B - D))
if dists:
    print(f'  B-only variant distance from TSS: min {dists[0]:,}  median {dists[len(dists) // 2]:,}  max {dists[-1]:,}')

def stat(files, want):
    for fp in files:
        with open(fp) as fh:
            h = fh.readline().rstrip('\n').split('\t')
            if 'Pvalue' not in h:
                continue
            idx = {c: h.index(c) for c in h}
            for line in fh:
                fl = line.rstrip('\n').split('\t')
                if fl[idx['Region']].split('.')[0] == gid.split('.')[0]:
                    return {c: fl[idx[c]] for c in ('Pvalue', 'Pvalue_Burden', 'Pvalue_SKAT', 'BETA_Burden', 'SE_Burden', 'MAC', 'Number_rare') if c in idx}
    return None
print('\n=== AP1M2 test statistics ===')
for arm, fs in (('B', sorted(glob.glob(f'{R}/saige_step2_pilot/tchl.chr19.B_3kb_re2g.B.part[0-9][0-9][0-9]'))), ('D', [f'{R}/saige_step2_pilot/tchl.chr19.D_3kb_only']), ('C', [f'{R}/saige_step2_pilot/tchl.chr19.C_nearest']), ('A', sorted(glob.glob(f'{R}/saige_step2_v4/tchl.chr19.part[0-9][0-9][0-9]')))):
    s = stat(fs, gid)
    print(f'  {arm}: {s}')
links = 0
tissues = collections.Counter()
for fp in sorted(glob.glob(f'{R}/re2g/*.bed.gz')):
    with gzip.open(fp, 'rt') as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 15:
                continue
            if f[13].split('.')[0] == gid.split('.')[0]:
                links += 1
                tissues[fp.split('/')[-1]] += 1
print(f'\nrE2G links targeting AP1M2: {links} across {len(tissues)} tissue files')
for k, v in tissues.most_common(5):
    print(f'   {k}: {v}')
