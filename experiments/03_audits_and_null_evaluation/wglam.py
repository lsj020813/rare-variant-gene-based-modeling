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
    for _ in range(90):
        mid = (lo + hi) / 2
        if erfc(sqrt(mid / 2)) > p:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2
O = _configured('${PROJECT_ROOT}/work/ref/saige_step2_v4')
THR = 2.5e-06
print(f"{'trait':<6}{'genes':>8}{'lambda':>9}{'min_p':>11}{'sig':>6}{'p<1e-4':>8}")
res = {}
for t in ('htn', 'dm', 'lip', 'tchl'):
    ps = []
    sig = []
    for f in sorted(glob.glob(f'{O}/{t}.chr*.part[0-9][0-9][0-9]')):
        with open(f) as fh:
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
                if 0 < p <= 1:
                    ps.append(p)
                    if p < THR:
                        sig.append((fl[gi], p))
    lam = median([qchisq1(p) for p in ps]) / 0.4549364
    res[t] = (len(ps), lam, min(ps), sig)
    print(f'{t:<6}{len(ps):>8,}{lam:>9.3f}{min(ps):>11.2e}{len(sig):>6}{sum((1 for p in ps if p < 0.0001)):>8}')
print()
for t, (n, lam, mn, sig) in res.items():
    if sig:
        print(f'{t} significant genes ({len(sig)}):')
        for g, p in sorted(sig, key=lambda x: x[1])[:12]:
            print(f'   {g:<22} {p:.2e}')
