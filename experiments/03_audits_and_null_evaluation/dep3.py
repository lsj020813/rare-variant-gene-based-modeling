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
import gzip, subprocess, json
R = _configured('${PROJECT_ROOT}/work/ref')
BCF = _configured('${BCFTOOLS}')

def norm(k):
    return k[3:] if k.startswith('chr') else k
sig = json.load(open(f'{R}/l1/tchl_B_full_verdict.json'))['sig']
span = {}
sym = {}
with gzip.open(f'{R}/deductive/gencode.sorted.gtf.gz', 'rt') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.split('\t', 9)
        if len(f) < 9 or f[2] != 'gene':
            continue
        a = f[8]
        gid = a.split('gene_id "')[1].split('"')[0].split('.')[0]
        span[gid] = (f[0], int(f[3]) - 3000, int(f[4]) + 3000)
        if 'gene_name "' in a:
            sym[gid] = a.split('gene_name "')[1].split('"')[0]
posmap = {}

def load_chrom(chn):
    if chn in posmap:
        return posmap[chn]
    m = {}
    out = subprocess.run(f"{BCF} query -f '%POS\t%ID\n' {R}/lift38_keyed/chr{chn}.keyed38.vcf.gz", shell=True, capture_output=True, text=True)
    for l in out.stdout.splitlines():
        p, kid = l.split('\t')
        m[norm(kid)] = int(p)
    posmap[chn] = m
    return m
dep = []
for g, p in sorted(sig.items(), key=lambda x: x[1]):
    if g not in span:
        continue
    ch38, lo, hi = span[g]
    chn = ch38.replace('chr', '')
    keys = []
    for line in open(f'{R}/groupfiles_bwg/chr{chn}.B_3kb_re2g.txt'):
        f = line.split()
        if f[0].split('.')[0] == g and f[1] == 'var':
            keys = [norm(k) for k in f[2:]]
            break
    m = load_chrom(chn)
    mapped = {k: m[k] for k in keys if k in m}
    assert len(mapped) >= 0.95 * len(keys), f'GATE FAIL {g}: mapped {len(mapped)}/{len(keys)} — key format mismatch'
    n_dist = sum((1 for pos in mapped.values() if lo <= pos <= hi))
    n_link = len(mapped) - n_dist
    dep.append((sym.get(g, g), p, len(keys), len(mapped), n_dist, n_link))
    print(f'{sym.get(g, g):<10} p={p:.1e}  vars {len(keys):>4} mapped {len(mapped):>4} = {n_dist:>4} dist + {n_link:>4} link ({n_link / max(len(mapped), 1) * 100:.0f}% link)', flush=True)
json.dump(dep, open(f'{R}/l1/tchl_dep_corrected.json', 'w'), indent=1)
print('DEP2_DONE')
