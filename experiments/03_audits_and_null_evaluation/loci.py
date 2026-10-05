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
import glob, gzip, json, collections
R = _configured('${PROJECT_ROOT}/work/ref')
THR = 2.5e-06
coord = {}
with gzip.open(f'{R}/deductive/gencode.sorted.gtf.gz', 'rt') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.split('\t')
        if len(f) < 9 or f[2] != 'gene':
            continue
        if f[0] not in ('19', 'chr19'):
            continue
        a = f[8]
        if 'gene_type "protein_coding"' not in a:
            continue
        gid = a.split('gene_id "')[1].split('"')[0]
        coord[gid.split('.')[0]] = (int(f[3]), int(f[4]), f[6], a.split('gene_name "')[1].split('"')[0] if 'gene_name "' in a else '?')
print(f'chr19 protein-coding genes in GTF: {len(coord):,}')

def read_res(files):
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
                g = fl[gi]
                if g not in out or p < out[g]:
                    out[g] = p
    return out

def read_group(files):
    out = collections.defaultdict(set)
    for fp in files:
        for line in open(fp):
            f = line.split()
            if len(f) > 2 and f[1] == 'var':
                out[f[0]].update(f[2:])
    return out
ARMS = {'A_100kb_re2g': (sorted(glob.glob(f'{R}/saige_step2_v4/tchl.chr19.part[0-9][0-9][0-9]')), sorted(glob.glob(f'{R}/groupfiles_chunks/chr19.part*.txt'))), 'B_3kb_re2g': (sorted(glob.glob(f'{R}/saige_step2_pilot/tchl.chr19.B_3kb_re2g.B.part[0-9][0-9][0-9]')), [f'{R}/groupfiles_pilot/chr19.B_3kb_re2g.txt']), 'C_nearest': ([f'{R}/saige_step2_pilot/tchl.chr19.C_nearest'], [f'{R}/groupfiles_pilot/chr19.C_nearest.txt']), 'D_3kb_only': ([f'{R}/saige_step2_pilot/tchl.chr19.D_3kb_only'], [f'{R}/groupfiles_pilot/chr19.D_3kb_only.txt'])}

def cluster_dist(sig, win=500000):
    pts = []
    for g, p in sig.items():
        b = g.split('.')[0]
        if b in coord:
            pts.append((coord[b][0], g, p))
    pts.sort()
    loci = []
    cur = []
    for pos, g, p in pts:
        if cur and pos - cur[-1][0] > win:
            loci.append(cur)
            cur = []
        cur.append((pos, g, p))
    if cur:
        loci.append(cur)
    return (loci, len(pts))

def cluster_share(sig, gv, frac=0.5):
    gs = list(sig)
    parent = {g: g for g in gs}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for i in range(len(gs)):
        for j in range(i + 1, len(gs)):
            a, b = (gv.get(gs[i], set()), gv.get(gs[j], set()))
            if not a or not b:
                continue
            if len(a & b) / min(len(a), len(b)) >= frac:
                ra, rb = (find(gs[i]), find(gs[j]))
                if ra != rb:
                    parent[ra] = rb
    grp = collections.defaultdict(list)
    for g in gs:
        grp[find(g)].append(g)
    return list(grp.values())
report = {}
for arm, (rfs, gfs) in ARMS.items():
    res = read_res(rfs)
    gv = read_group(gfs)
    sig = {g: p for g, p in res.items() if p < THR}
    ld, _ = cluster_dist(sig)
    ls = cluster_share(sig, gv)
    leads = []
    for L in ld:
        best = min(L, key=lambda x: x[2])
        leads.append((best[1], best[2], len(L)))
    report[arm] = dict(tested=len(res), sig=len(sig), loci_dist=len(ld), loci_share=len(ls), spread_dist=round(len(sig) / len(ld), 2) if ld else 0, spread_share=round(len(sig) / len(ls), 2) if ls else 0, leads=[(g, p, n) for g, p, n in sorted(leads, key=lambda x: x[1])])
    print(f'\n=== {arm} ===')
    print(f'  tested {len(res):,}  significant {len(sig)}')
    print(f"  loci (거리 500kb): {len(ld)}   좌위당 유전자 {report[arm]['spread_dist']}")
    print(f"  loci (변이 50% 공유): {len(ls)}   좌위당 유전자 {report[arm]['spread_share']}")
    for g, p, n in sorted(leads, key=lambda x: x[1])[:6]:
        nm = coord.get(g.split('.')[0], ('', '', '', '?'))[3]
        print(f'    lead {nm:<12} {g:<20} p={p:.2e}  ({n} genes in locus)')
json.dump(report, open(_configured('${PROJECT_ROOT}/work/run_audit/locus_report.json'), 'w'), indent=1)
print('\nwritten locus_report.json')
