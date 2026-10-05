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
import glob
from statistics import median
from math import sqrt, erfc

def qchisq1(p):
    lo, hi = (0.0, 1000.0)
    for _ in range(80):
        mid = (lo + hi) / 2
        if erfc(sqrt(mid / 2)) > p:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2

def summarize(files, tag):
    ps = []
    n = 0
    for fp in files:
        with open(fp) as fh:
            h = fh.readline().rstrip('\n').split('\t')
            if 'Pvalue' not in h:
                continue
            pi = h.index('Pvalue')
            for line in fh:
                f = line.rstrip('\n').split('\t')
                n += 1
                try:
                    p = float(f[pi])
                except (ValueError, IndexError):
                    continue
                if 0 < p <= 1:
                    ps.append(p)
    if not ps:
        print(f'  {tag}: none')
        return
    lam = median([qchisq1(p) for p in ps]) / 0.4549364
    print(f'  {tag:<26} n={n:<5} lambda={lam:.3f} min_p={min(ps):.2e} p<0.05={sum((1 for p in ps if p < 0.05))}  p<1e-4={sum((1 for p in ps if p < 0.0001))}')
print('=== v4 vs v2 — chr22 ===')
for t in ('htn', 'dm', 'lip', 'tchl'):
    summarize(sorted(glob.glob(_configured(f'${{PROJECT_ROOT}}/work/ref/saige_step2_v4/{t}.chr22.part[0-9][0-9][0-9]'))), f'v4 {t}')
    summarize(sorted(glob.glob(_configured(f'${{PROJECT_ROOT}}/work/ref/saige_step2_v2/{t}.chr22.part[0-9][0-9][0-9]'))), f'v2 {t}')
