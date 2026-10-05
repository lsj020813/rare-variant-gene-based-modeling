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
import glob, gzip, json, collections, re
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
part_re = re.compile('\\.part\\d{3}$')
res = {}
for T in ('tchl', 'htn', 'dm', 'lip'):
    sig = {}
    for fp in glob.glob(f'{R}/saige_step2_v4/{T}.chr*.part*'):
        if not part_re.search(fp):
            continue
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
                if p < THR:
                    sig[f[gi].split('.')[0]] = p
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
    res[T] = {'n_sig': len(sig), 'n_loci': len(locs), 'loci': [{'chrom': L['chrom'], 'top': sym.get(L['top'], L['top']), 'topp': L['topp'], 'n_genes': len(L['genes'])} for L in locs]}
    print(f'{T}: sig {len(sig)}  loci {len(locs)}', flush=True)
json.dump(res, open(f'{R}/l1/v4_verdict_side.json', 'w'), indent=1)
print('V4_SIDE_DONE')
