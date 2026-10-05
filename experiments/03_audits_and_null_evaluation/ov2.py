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
G = _configured('${PROJECT_ROOT}/work/ref/groupfiles_chunks')
GENES = {'ENSG00000130203': 'APOE', 'ENSG00000130204': 'TOMM40', 'ENSG00000130208': 'APOC1'}
sets = {}
for fp in glob.glob(f'{G}/chr19.part*.txt'):
    for line in open(fp):
        f = line.split()
        if len(f) > 2 and f[1] == 'var':
            g = f[0].split('.')[0]
            if g in GENES:
                sets[GENES[g]] = set(f[2:])
print('=== 겹침은 그룹파일의 성질인가? (널모델과 무관) ===')
ks = list(sets)
for i in range(len(ks)):
    for j in range(i + 1, len(ks)):
        a, b = (sets[ks[i]], sets[ks[j]])
        if not a or not b:
            continue
        print(f'  {ks[i]:<8} vs {ks[j]:<8} shared={len(a & b):<5} jaccard={len(a & b) / len(a | b):.3f}  (|A|={len(a)} |B|={len(b)})')
print()
print('=== v4 에서 chr19 가 끝났나? ===')
done = sorted(glob.glob(_configured('${PROJECT_ROOT}/work/ref/saige_step2_v4/*.chr19.part*.done')))
print(f"  chr19 done chunks: {len(done)} / {len(glob.glob(f'{G}/chr19.part*.txt')) * 4}")
