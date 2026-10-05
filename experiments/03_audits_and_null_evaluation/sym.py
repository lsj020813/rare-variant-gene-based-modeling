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
THR = 2.5e-06
sym = {}
with gzip.open(f'{R}/deductive/gencode.sorted.gtf.gz', 'rt') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.split('\t', 9)
        if len(f) < 9 or f[2] != 'gene':
            continue
        a = f[8]
        if 'gene_name "' not in a or 'gene_id "' not in a:
            continue
        sym[a.split('gene_id "')[1].split('"')[0].split('.')[0]] = (a.split('gene_name "')[1].split('"')[0], f[0], int(f[3]))
print(f'GTF genes with symbols: {len(sym):,}')
O = f'{R}/saige_step2_v4'
for t in ('tchl', 'htn', 'dm', 'lip'):
    hits = []
    for fp in sorted(glob.glob(f'{O}/{t}.chr*.part[0-9][0-9][0-9]')):
        with open(fp) as fh:
            h = fh.readline().rstrip('\n').split('\t')
            if 'Pvalue' not in h:
                continue
            gi = h.index('Region')
            pi = h.index('Pvalue')
            ni = h.index('Number_rare')
            bu = h.index('Pvalue_Burden')
            be = h.index('BETA_Burden')
            for line in fh:
                fl = line.rstrip('\n').split('\t')
                try:
                    p = float(fl[pi])
                except (ValueError, IndexError):
                    continue
                if p < THR:
                    g = fl[gi].split('.')[0]
                    nm, ch, pos = sym.get(g, ('?', '?', 0))
                    hits.append((p, nm, ch, pos, fl[ni], fl[bu], fl[be]))
    hits.sort()
    print(f'\n===== {t.upper()} — {len(hits)} significant genes =====')
    for p, nm, ch, pos, n, bp, bb in hits:
        try:
            bpf = f'{float(bp):.3f}'
        except ValueError:
            bpf = bp
        print(f"  {nm:<14} chr{ch.replace('chr', '')}:{pos / 1000000.0:>7.2f}Mb  p={p:.2e}  n={n:<5} burden_p={bpf}")
