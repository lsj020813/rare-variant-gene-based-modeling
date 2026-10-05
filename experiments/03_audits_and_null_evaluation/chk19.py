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
R = _configured('${PROJECT_ROOT}/work/ref')

def read(fp_pat):
    out = {}
    for fp in sorted(glob.glob(fp_pat)):
        with open(fp) as fh:
            h = fh.readline().rstrip('\n').split('\t')
            if 'Pvalue' not in h:
                continue
            gi = h.index('Region')
            pi = h.index('Pvalue')
            bu = h.index('Pvalue_Burden')
            sk = h.index('Pvalue_SKAT')
            be = h.index('BETA_Burden')
            se = h.index('SE_Burden')
            mc = h.index('MAC')
            for line in fh:
                f = line.rstrip('\n').split('\t')
                try:
                    out[f[gi]] = (f[pi], f[bu], f[sk], f[be], f[se], f[mc])
                except IndexError:
                    pass
    return out
pilot = read(f'{R}/saige_step2_pilot/tchl.chr19.B_3kb_re2g.B.part[0-9][0-9][0-9]')
wg = read(f'{R}/saige_step2_bwg/tchl.chr19')
print(f'pilot genes: {len(pilot):,}   wg genes: {len(wg):,}')
common = set(pilot) & set(wg)
only_p = set(pilot) - set(wg)
only_w = set(wg) - set(pilot)
print(f'common: {len(common):,}  pilot-only: {len(only_p)}  wg-only: {len(only_w)}')
mm = []
for g in common:
    for i, cn in enumerate(('Pvalue', 'Burden', 'SKAT', 'BETA', 'SE', 'MAC')):
        if pilot[g][i] != wg[g][i]:
            mm.append((g, cn, pilot[g][i], wg[g][i]))
print(f'field mismatches: {len(mm)} / {len(common) * 6:,}')
for m in mm[:8]:
    print('  ', m)
print('VERDICT:', 'IDENTICAL' if not mm and (not only_p) and (not only_w) else 'MISMATCH')
