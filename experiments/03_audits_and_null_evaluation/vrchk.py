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
import csv, subprocess, os
P = _configured('${ANALYSIS_ROOT}/work/saige_gene/varratio_plink/${VARIANCE_RATIO_PLINK_BASENAME}')
Q = _configured('${ANALYSIS_ROOT}/work/saige_gene/varratio_plink_chr22_gp90/chr22.GP90.varratio_mac10plus_retry3')
for tag, p in [('mac10_2000', P), ('gp90_retry3', Q)]:
    n_mark = sum((1 for _ in open(p + '.bim')))
    n_samp = sum((1 for _ in open(p + '.fam')))
    print(f'[{tag}] markers {n_mark:,} samples {n_samp:,}')
for tag, p in [('mac10_2000', P)]:
    ac = p + '.freq_counts.acount'
    if os.path.exists(ac):
        lo = hi = 0
        with open(ac) as fh:
            h = fh.readline().split()
            try:
                i = h.index('ALT_CTS')
            except ValueError:
                i = 4
            for line in fh:
                f = line.split()
                try:
                    v = float(f[i])
                except (ValueError, IndexError):
                    continue
                if 10 < v <= 20.5:
                    lo += 1
                elif v > 20.5:
                    hi += 1
        print(f'  [{tag}] full-cohort MAC 10-20.5: {lo:,} | >20.5: {hi:,}')
    else:
        print(f'  [{tag}] no acount file at {ac}')
ph = {r.split('\t')[0] for i, r in enumerate(open(_configured('${PROJECT_ROOT}/work/ref/pheno_v3/htn_v3.tsv'))) if i}
for tag, p in [('mac10_2000', P), ('gp90_retry3', Q)]:
    fam = {l.split()[1] for l in open(p + '.fam')}
    print(f'  [{tag}] phenotype samples -> in fileset: {len(ph & fam):,}')
