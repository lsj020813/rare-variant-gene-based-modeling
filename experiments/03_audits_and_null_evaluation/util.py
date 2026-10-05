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
import csv, subprocess, glob, os
BCF = _configured('${BCFTOOLS}')
gt = {s.split('_')[0] for s in subprocess.run(_configured(f'{BCF} query -l ${{PROJECT_ROOT}}/work/ref/band_vcf/chr22.band.vcf.gz'), shell=True, capture_output=True, text=True).stdout.split()}
print(f'genotyped (denominator): {len(gt):,}')

def tab(p):
    with open(p, encoding='utf-8', errors='ignore') as fh:
        rd = csv.reader(fh, delimiter='\t')
        hdr = next(rd)
        rows = [r for r in rd if r]
    return (hdr, rows)
P = _configured('${PROJECT_ROOT}/work/ref/pheno_v2')
T = _configured('${ANALYSIS_ROOT}/work/phenotype/tchl_primary_complete_case.tsv')
print()
print(f"{'trait':<7}{'in table':>10}{'analysed':>10}{'cases':>9}{'% of gt':>9}")
for name, path, ycol in [('HTN', f'{P}/htn_v2.tsv', 'y'), ('DM', f'{P}/dm_v2.tsv', 'y'), ('LIP', f'{P}/lip_v2.tsv', 'y'), ('TCHL', T, 'TCHL_rint')]:
    hdr, rows = tab(path)
    yi = hdr.index(ycol)
    ids = {r[0].split('_')[0] for r in rows}
    ok = 0
    cases = 0
    for r in rows:
        if r[0].split('_')[0] not in gt:
            continue
        v = r[yi].strip()
        if v in ('', 'NA', 'NaN', '.'):
            continue
        ok += 1
        try:
            if float(v) == 1.0:
                cases += 1
        except ValueError:
            pass
    print(f"{name:<7}{len(ids):>10,}{ok:>10,}{(cases if name != 'TCHL' else 0):>9,}{ok / len(gt) * 100:>8.1f}%")
