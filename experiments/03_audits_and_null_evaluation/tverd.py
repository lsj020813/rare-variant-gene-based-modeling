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
import glob, gzip, json, collections, os
R = _configured('${PROJECT_ROOT}/work/ref')
THR = 2.5e-06
coord = {}
sym = {}
with gzip.open(f'{R}/deductive/gencode.sorted.gtf.gz', 'rt') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.split('\t', 9)
        if len(f) < 9 or f[2] != 'gene':
            continue
        a = f[8]
        if 'gene_type "protein_coding"' not in a:
            continue
        gid = a.split('gene_id "')[1].split('"')[0].split('.')[0]
        coord[gid] = (f[0], int(f[3]))
        if 'gene_name "' in a:
            sym[gid] = a.split('gene_name "')[1].split('"')[0]

def read_sig(paths):
    sig = {}
    for fp in paths:
        with open(fp) as fh:
            h = fh.readline().rstrip('\n').split('\t')
            if 'Pvalue' not in h:
                continue
            gi, pi = (h.index('Region'), h.index('Pvalue'))
            for line in fh:
                f = line.rstrip('\n').split('\t')
                try:
                    p = float(f[pi])
                except (ValueError, IndexError):
                    continue
                g = f[gi].split('.')[0]
                if p < THR and (g not in sig or p < sig[g]):
                    sig[g] = p
    return sig
paths = [f'{R}/saige_step2_bwg/tchl.chr{i}' for i in range(2, 23) if os.path.getsize(f'{R}/saige_step2_bwg/tchl.chr{i}.done') > 0]
paths += sorted(glob.glob(f'{R}/saige_step2_bwg_chr1/tchl.chr1.part*[0-9]'))
paths = [p for p in paths if not p.endswith(('.log', '.done'))]
sig = read_sig(paths)

def cluster(sig):
    locs = []
    for g, p in sorted(sig.items(), key=lambda x: x[1]):
        if g not in coord:
            continue
        ch, ts = coord[g]
        for L in locs:
            if L['chrom'] == ch and abs(ts - L['tss']) < 500000:
                L['genes'].append(g)
                break
        else:
            locs.append({'chrom': ch, 'tss': ts, 'genes': [g], 'top': g, 'topp': p})
    return locs
locs = cluster(sig)
print(f'B arm TCHL whole genome: sig {len(sig)}  loci {len(locs)}')
for L in locs:
    print(f"  {L['chrom']:>6} {sym.get(L['top'], L['top']):<10} p={L['topp']:.1e} genes={len(L['genes'])}")
gene_span = {}
with gzip.open(f'{R}/deductive/gencode.sorted.gtf.gz', 'rt') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.split('\t', 9)
        if len(f) < 9 or f[2] != 'gene':
            continue
        a = f[8]
        gid = a.split('gene_id "')[1].split('"')[0].split('.')[0]
        gene_span[gid] = (f[0].replace('chr', ''), int(f[3]) - 3000, int(f[4]) + 3000)
dep = []
for g, p in sorted(sig.items(), key=lambda x: x[1]):
    ch = coord[g][0].replace('chr', '') if g in coord else None
    gf = f'{R}/groupfiles_bwg/chr{ch}.B_3kb_re2g.txt'
    if not os.path.exists(gf):
        continue
    vars_g = []
    for line in open(gf):
        f = line.split()
        if f[0].split('.')[0] == g and f[1] == 'var':
            vars_g = f[2:]
            break
    if not vars_g or g not in gene_span:
        continue
    _, lo, hi = gene_span[g]
    n_dist = sum((1 for v in vars_g if lo <= int(v.split(':')[1]) <= hi))
    n_link = len(vars_g) - n_dist
    dep.append((sym.get(g, g), p, len(vars_g), n_dist, n_link))
print('\nsig genes: variants total / within-3kb / rE2G-only')
for s, p, nt, nd, nl in dep:
    print(f'  {s:<10} p={p:.1e}  {nt:>4} = {nd:>4} dist + {nl:>4} link ({nl / nt * 100:.0f}% link)')
json.dump({'sig': {k: v for k, v in sig.items()}, 'loci': [{'chrom': L['chrom'], 'top': sym.get(L['top'], L['top']), 'topp': L['topp'], 'n': len(L['genes'])} for L in locs], 'dep': dep}, open(f'{R}/l1/tchl_B_full_verdict.json', 'w'), indent=1)
print('TCHL_VERDICT_DONE')
