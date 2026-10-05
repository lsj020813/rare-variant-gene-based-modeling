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
import glob, os
O = _configured('${PROJECT_ROOT}/work/ref/saige_step2_v4')
G = _configured('${PROJECT_ROOT}/work/ref/groupfiles_chunks')
traits = ('htn', 'dm', 'lip', 'tchl')
miss = []
short = []
ok = 0
tot_rows = {t: 0 for t in traits}
for ch in range(1, 23):
    for gf in sorted(glob.glob(f'{G}/chr{ch}.part*.txt')):
        part = os.path.basename(gf)[:-4].split('.', 1)[1]
        exp = len({l.split()[0] for l in open(gf) if len(l.split()) > 2 and l.split()[1] == 'var'})
        for t in traits:
            f = f'{O}/{t}.chr{ch}.{part}'
            if not os.path.exists(f) or os.path.getsize(f) == 0:
                miss.append(f'{t}.chr{ch}.{part}')
                continue
            n = sum((1 for _ in open(f))) - 1
            tot_rows[t] += n
            if n != exp:
                short.append(f'{t}.chr{ch}.{part} {n}/{exp}')
            else:
                ok += 1
print(f"expected runs: {len(glob.glob(f'{G}/chr*.part*.txt')) * 4}")
print(f'  exact-row OK : {ok}')
print(f'  missing/empty: {len(miss)}')
print(f'  row mismatch : {len(short)}')
for m in (miss + short)[:6]:
    print('   ', m)
print('\ngene-tests per trait:')
for t in traits:
    print(f'  {t:<5} {tot_rows[t]:,}')
