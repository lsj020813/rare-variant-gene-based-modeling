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
import glob, gzip
R = _configured('${PROJECT_ROOT}/work/ref')
RV = ['PCSK9', 'APOB', 'ANGPTL3', 'ANGPTL4', 'APOC3', 'LPL', 'ABCA1', 'LCAT', 'LIPG', 'SCARB1', 'LDLR', 'APOE', 'CETP', 'NPC1L1', 'MYLIP', 'STAP1', 'LIPC', 'SREBF1']
CV = ['APOE', 'LDLR', 'CETP', 'HMGCR', 'SORT1', 'CELSR2', 'PSRC1', 'APOA1', 'APOA5', 'TRIB1', 'ABCG8', 'LPA', 'TIMD4', 'MLXIPL', 'GCKR', 'DOCK7', 'ANGPTL8']
want = sorted(set(RV + CV))
sym = {}
with gzip.open(f'{R}/deductive/gencode.sorted.gtf.gz', 'rt') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.split('\t', 9)
        if len(f) < 9 or f[2] != 'gene':
            continue
        a = f[8]
        if 'gene_name "' not in a:
            continue
        nm = a.split('gene_name "')[1].split('"')[0]
        if nm in want:
            sym[a.split('gene_id "')[1].split('"')[0].split('.')[0]] = (nm, f[0])
res = {}
for fp in sorted(glob.glob(f'{R}/saige_step2_v4/tchl.chr*.part[0-9][0-9][0-9]')):
    with open(fp) as fh:
        h = fh.readline().rstrip('\n').split('\t')
        if 'Pvalue' not in h:
            continue
        gi = h.index('Region')
        pi = h.index('Pvalue')
        ni = h.index('Number_rare')
        bu = h.index('Pvalue_Burden')
        for line in fh:
            fl = line.rstrip('\n').split('\t')
            g = fl[gi].split('.')[0]
            if g in sym:
                try:
                    res[sym[g][0]] = (float(fl[pi]), fl[ni], fl[bu], sym[g][1])
                except (ValueError, IndexError):
                    pass
THR = 2.5e-06

def show(title, names):
    print(f'\n=== {title} ===')
    hit = 0
    for nm in names:
        if nm in res:
            p, n, bp, ch = res[nm]
            mk = ' ***SIG' if p < THR else ''
            if p < THR:
                hit += 1
            print(f"  {nm:<10} chr{ch.replace('chr', ''):<3} p={p:.2e}  n={n:<5} burden={float(bp):.3f}{mk}")
        else:
            print(f'  {nm:<10} NOT TESTED')
    print(f'  -> significant: {hit}/{len(names)}')
show('RARE-variant (sequencing) lipid genes', sorted(set(RV)))
show('COMMON-variant (array GWAS) lipid loci', sorted(set(CV)))
