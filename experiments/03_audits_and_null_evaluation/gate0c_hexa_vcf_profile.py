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
import os, re, subprocess, collections, json
BCF = _configured('${BCFTOOLS}')
ROOT = _configured('${COHORT_DATA}')
V = _required('HEXA_VCF_FILE')

def shape(s):
    return re.sub('[0-9]', 'N', s)
out = {'vcf_size_gb': round(os.path.getsize(V) / 1000000000.0, 2), 'index_present': [x for x in ('.tbi', '.csi') if os.path.exists(V + x)]}
hdr = subprocess.run([BCF, 'view', '-h', V], capture_output=True, text=True, check=True).stdout.splitlines()
out['n_header_lines'] = len(hdr)
out['INFO_fields'] = [re.search('ID=([^,]+)', l).group(1) for l in hdr if l.startswith('##INFO')]
out['FORMAT_fields'] = [re.search('ID=([^,]+)', l).group(1) for l in hdr if l.startswith('##FORMAT')]
out['contig_lines'] = sum((1 for l in hdr if l.startswith('##contig')))
out['ref_or_build_lines'] = [l[:120] for l in hdr if re.search('reference|assembly|GRCh|hg19|hg38|source|imput|minimac|beagle|eagle|panel', l, re.I)][:8]
ids = set(subprocess.run([BCF, 'query', '-l', V], capture_output=True, text=True, check=True).stdout.split())
out['vcf_n_samples'] = len(ids)
out['id_shapes'] = collections.Counter((shape(x) for x in ids)).most_common(3)
for k, d in {'HEXA_850K': _required('HEXA_IDAT_DIR'), 'KARE_850K': _required('KARE_HIGH_DENSITY_IDAT_DIR'), 'KARE_450K': _required('KARE_LOW_DENSITY_IDAT_DIR')}.items():
    stems = set((f.rsplit('_', 1)[0] for f in os.listdir(os.path.join(ROOT, d)) if f.endswith('.idat')))
    out[k + '_in_hexa_vcf'] = len(stems & ids)
try:
    out['chr22_records_index'] = subprocess.run([BCF, 'index', '-n', V], capture_output=True, text=True).stdout.strip()
except Exception as e:
    out['chr22_records_index'] = 'NA'
print(json.dumps(out, ensure_ascii=False, indent=1))
