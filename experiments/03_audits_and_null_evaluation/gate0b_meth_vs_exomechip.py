#!/usr/bin/env python
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
import os, re, collections, json
ROOT = _configured('${COHORT_DATA}')
ped = _required('EXOME_PED_FILE')
fid = set()
iid = set()
with open(ped) as f:
    for line in f:
        p = line.split(None, 2)
        fid.add(p[0])
        iid.add(p[1])

def shape(s):
    return re.sub('[0-9]', 'N', s)
out = {'ped_n': len(iid), 'iid_shapes': collections.Counter((shape(x) for x in iid)).most_common(3), 'fid_eq_iid': len(fid & iid)}
for k, d in {'KARE_850K': _required('KARE_HIGH_DENSITY_IDAT_DIR'), 'KARE_450K': _required('KARE_LOW_DENSITY_IDAT_DIR'), 'HEXA_850K': _required('HEXA_IDAT_DIR')}.items():
    stems = set((f.rsplit('_', 1)[0] for f in os.listdir(os.path.join(ROOT, d)) if f.endswith('.idat')))
    out[k + '_in_exomechip_iid'] = len(stems & iid)
k8 = set((f.rsplit('_', 1)[0] for f in os.listdir(_required('KARE_HIGH_DENSITY_IDAT_DIR')) if f.endswith('.idat')))
k4 = set((f.rsplit('_', 1)[0] for f in os.listdir(_required('KARE_LOW_DENSITY_IDAT_DIR')) if f.endswith('.idat')))
out['KARE_850K_and_450K_overlap'] = len(k8 & k4)
csv = _required('METHYLATION_FILE')
with open(csv) as f:
    hdr = f.readline().rstrip('\n').split(',')
out['csv50_ncols'] = len(hdr)
out['csv50_hdr_shapes'] = collections.Counter((shape(h) for h in hdr)).most_common(3)
print(json.dumps(out, ensure_ascii=False, indent=1))
