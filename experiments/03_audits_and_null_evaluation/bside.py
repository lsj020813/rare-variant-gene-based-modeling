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
res = {}
for T in ('tchl', 'htn', 'dm', 'lip'):
    sig = {}
    lam_ps = []
    for ch in list(range(2, 23)):
        fp = f'{R}/saige_step2_bwg/{T}.chr{ch}'
        if not (os.path.exists(fp + '.done') and os.path.getsize(fp + '.done') > 0):
            continue
        with open(fp) as fh:
            h = fh.readline().rstrip('\n').split('\t')
            gi, pi = (h.index('Region'), h.index('Pvalue'))
            for line in fh:
                f = line.rstrip('\n').split('\t')
                try:
                    p = float(f[pi])
                except (ValueError, IndexError):
                    continue
                lam_ps.append(p)
                if p < THR:
                    sig[f[gi].split('.')[0]] = p
    from math import sqrt, erfc

    def qc(p):
        lo, hi = (0.0, 1000.0)
        for _ in range(80):
            mid = (lo + hi) / 2
            if erfc(sqrt(mid / 2)) > p:
                lo = mid
            else:
                hi = mid
        return lo
    lam_ps.sort()
    med = lam_ps[len(lam_ps) // 2]
    lam = qc(med) / 0.4549
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
    res[T] = {'n_tests': len(lam_ps), 'lambda': round(lam, 3), 'n_sig': len(sig), 'n_loci': len(locs), 'loci': [{'chrom': L['chrom'], 'top': sym.get(L['top'], L['top']), 'topp': L['topp'], 'n_genes': len(L['genes'])} for L in locs]}
    print(f'{T}: tests {len(lam_ps):,} lambda {lam:.3f} sig {len(sig)} loci {len(locs)}', flush=True)
    for L in res[T]['loci']:
        print('   ', L['chrom'], L['top'], f"p={L['topp']:.1e}", f"genes={L['n_genes']}")
json.dump(res, open(f'{R}/l1/bwg_verdict_side.json', 'w'), indent=1)
print('BWG_SIDE_DONE (chr2-22; chr1 pending)')
