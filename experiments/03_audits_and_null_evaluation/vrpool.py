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
import subprocess, os
P2 = _configured('${PLINK2}')
sets = {'gp90_retry3': _configured('${PROJECT_ROOT}/work/ref/vr_src/gp90'), 'mac10_2000': _configured('${PROJECT_ROOT}/work/ref/vr_src/mac2000')}
for tag, p in [('mac10_2000', _configured('${PROJECT_ROOT}/work/ref/annotation_ref'))]:
    pass
import glob
print('=== surviving PLINK filesets on disk ===')
for pat in [_configured('${DATA_ROOT}/GWAS/99.old.Saige/sparseGRM/*.fam'), _configured('${PROJECT_ROOT}/work/ref/**/*.fam')]:
    for f in glob.glob(pat, recursive=True):
        n = sum((1 for _ in open(f)))
        m = f[:-4] + '.bim'
        nm = sum((1 for _ in open(m))) if os.path.exists(m) else 0
        print(f'  samples {n:>7,}  markers {nm:>9,}  {f}')
