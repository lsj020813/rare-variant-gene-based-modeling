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
import os, re, subprocess, collections, json, sys
BCF = _configured('${BCFTOOLS}')
VCF = _configured('${PROJECT_ROOT}/work/ref/band_vcf/chr22.band.vcf.gz')
ROOT = _configured('${COHORT_DATA}')
dirs = {'HEXA_850K': _required('HEXA_IDAT_DIR'), 'KARE_850K': _required('KARE_HIGH_DENSITY_IDAT_DIR'), 'KARE_450K': _required('KARE_LOW_DENSITY_IDAT_DIR')}

def shape(s):
    return re.sub('[0-9]', 'N', s)
vcf_ids = set(subprocess.run([BCF, 'query', '-l', VCF], capture_output=True, text=True, check=True).stdout.split())
out = {'vcf_n': len(vcf_ids), 'vcf_shapes': collections.Counter((shape(x) for x in vcf_ids)).most_common(3)}
hexa = set(open(_required('HEXA_ID_FILE')).read().split())
out['hexa_idtxt_n'] = len(hexa)
out['hexa_idtxt_in_vcf'] = len(hexa & vcf_ids)
allm = set()
for k, d in dirs.items():
    stems = set((f.rsplit('_', 1)[0] for f in os.listdir(os.path.join(ROOT, d)) if f.endswith('.idat')))
    grn = sum((1 for f in os.listdir(os.path.join(ROOT, d)) if f.endswith('_Grn.idat')))
    red = sum((1 for f in os.listdir(os.path.join(ROOT, d)) if f.endswith('_Red.idat')))
    out[k] = {'n_stems': len(stems), 'grn': grn, 'red': red, 'shapes': collections.Counter((shape(s) for s in stems)).most_common(2), 'in_vcf': len(stems & vcf_ids), 'in_hexa_idtxt': len(stems & hexa)}
    allm |= stems
out['meth_union_n'] = len(allm)
out['meth_union_in_vcf'] = len(allm & vcf_ids)
out['overlap_850K_pairs'] = len(set.intersection(*[set((f.rsplit('_', 1)[0] for f in os.listdir(os.path.join(ROOT, dirs[k])) if f.endswith('.idat'))) for k in ('HEXA_850K', 'KARE_850K')]))
print(json.dumps(out, ensure_ascii=False, indent=1))
