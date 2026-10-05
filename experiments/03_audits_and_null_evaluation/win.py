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
import glob, gzip, subprocess, statistics as st

def mstats(pat, tag):
    per = {}
    for fp in glob.glob(pat):
        for line in open(fp):
            f = line.split()
            if len(f) > 2 and f[1] == 'var':
                per[f[0]] = len(f) - 2
    if not per:
        print(f'  {tag:<34} (no files)')
        return
    v = sorted(per.values())
    print(f'  {tag:<34} genes={len(per):<6} median_m={st.median(v):<7.0f} p90={v[int(0.9 * len(v))]:<7} max={v[-1]}')
print('=== 지금 도는 전장 정답지 (groupfiles_chunks) ===')
mstats(_configured('${PROJECT_ROOT}/work/ref/groupfiles_chunks/chr19.part*.txt'), 'chr19 (running truth set)')
mstats(_configured('${PROJECT_ROOT}/work/ref/groupfiles_chunks/chr22.part*.txt'), 'chr22 (running truth set)')
print('=== 파일럿 후보 팔 (groupfiles_pilot) ===')
for arm in ('B_3kb_re2g', 'C_nearest', 'D_3kb_only'):
    mstats(_configured(f'${{PROJECT_ROOT}}/work/ref/groupfiles_pilot/chr19.{arm}.txt'), f'chr19 {arm}')
