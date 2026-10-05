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
import glob, statistics as st
G = _configured('${PROJECT_ROOT}/work/ref/groupfiles_chunks')
P = _configured('${PROJECT_ROOT}/work/ref/groupfiles_pilot')

def prof(pat, tag):
    per = {}
    for fp in glob.glob(pat):
        for line in open(fp):
            f = line.split()
            if len(f) > 2 and f[1] == 'var':
                per.setdefault(f[0], set()).update(f[2:])
    if not per:
        return None
    ms = sorted((len(v) for v in per.values()))
    allv = set().union(*per.values())
    total = sum(ms)
    return dict(tag=tag, genes=len(per), med=st.median(ms), uniq=len(allv), redundancy=total / len(allv))
rows = []
rows.append(prof(f'{G}/chr19.part*.txt', 'A  100kb u rE2G (running)'))
for arm, lbl in [('B_3kb_re2g', 'B  3kb u rE2G'), ('C_nearest', 'C  nearest gene'), ('D_3kb_only', 'D  3kb only')]:
    rows.append(prof(f'{P}/chr19.{arm}.txt', lbl))
print(f"{'arm':<26}{'genes':>7}{'median m':>10}{'variants':>10}{'redundancy':>12}")
for r in rows:
    if r:
        print(f"{r['tag']:<26}{r['genes']:>7}{r['med']:>10.0f}{r['uniq']:>10,}{r['redundancy']:>11.2f}x")
print()
print('=== 지금 런의 완료 상태 ===')
d = glob.glob(_configured('${PROJECT_ROOT}/work/ref/saige_step2_v4/*.done'))
import collections
ch = collections.Counter((f.split('.chr')[1].split('.')[0] for f in d))
print(f'  chunks {len(d)}/348   완료 청크 보유 염색체: {sorted(ch, key=lambda x: int(x))}')
