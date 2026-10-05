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
import glob, math
from statistics import median
O = _configured('${PROJECT_ROOT}/work/ref/saige_step2_v4')
G = _configured('${PROJECT_ROOT}/work/ref/groupfiles_chunks')
exp = set()
for fp in glob.glob(f'{G}/chr22.part*.txt'):
    for l in open(fp):
        f = l.split()
        if len(f) > 2 and f[1] == 'var':
            exp.add(f[0])
print(f'group-file genes: {len(exp)}')
THR = 2.5e-06
hdr_shown = False
for t in ('htn', 'dm', 'lip', 'tchl'):
    ps = []
    genes = set()
    cols = None
    for fp in sorted(glob.glob(f'{O}/{t}.chr22.part[0-9][0-9][0-9]')):
        with open(fp) as fh:
            h = fh.readline().rstrip('\n').split('\t')
            if cols is None:
                cols = h
            gi = h.index('Region')
            pi = h.index('Pvalue')
            mi = h.index('MAC') if 'MAC' in h else None
            for line in fh:
                f = line.rstrip('\n').split('\t')
                if len(f) <= max(gi, pi):
                    continue
                genes.add(f[gi])
                try:
                    p = float(f[pi])
                except ValueError:
                    continue
                if 0 < p <= 1:
                    ps.append(p)
    if not ps:
        print(f'{t}: NO P-VALUES')
        continue
    chi = [__import__('scipy.stats', fromlist=['chi2']).chi2.isf(p, 1) for p in ps] if False else None
    import statistics

    def qchisq1(p):
        from math import sqrt, log, erf
        lo, hi = (0.0, 1000.0)
        for _ in range(80):
            mid = (lo + hi) / 2
            from math import erfc
            if erfc(sqrt(mid / 2)) > p:
                lo = mid
            else:
                hi = mid
        return (lo + hi) / 2
    lam = median([qchisq1(p) for p in ps]) / 0.4549364
    sig = [p for p in ps if p < THR]
    print(f'{t:<5} genes={len(genes):<4} p-values={len(ps):<4} lambda={lam:.3f} min_p={min(ps):.2e} sig(<2.5e-6)={len(sig)}  missing_vs_groupfile={len(exp - genes)}')
    if not hdr_shown:
        print('  columns:', cols[:12])
        hdr_shown = True
