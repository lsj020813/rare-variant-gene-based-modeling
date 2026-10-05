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
import json, glob, os, numpy as np
W = _configured('${PROJECT_ROOT}/work')
d = json.load(open(f'{W}/run_meth/out/gate2_B86_vs_gwas.json'))
rows = [r for r in d['rows'] if not r['near_lead_any_trait']]
names = {}
for f in glob.glob(f'{W}/ref/**/*gene*', recursive=True)[:50]:
    if os.path.isfile(f) and os.path.getsize(f) < 200000000.0 and any((f.endswith(x) for x in ('.tsv', '.txt', '.bed', '.csv'))):
        try:
            with open(f) as fh:
                head = fh.readline()
                if 'ENSG' in head or True:
                    for line in fh:
                        if 'ENSG' in line:
                            q = line.replace(',', '\t').split('\t')
                            eid = [x for x in q if x.startswith('ENSG')]
                            sym = [x for x in q if x.isupper() and x.isalnum() and (not x.startswith('ENSG')) and (len(x) <= 12)]
                            if eid and sym:
                                names[eid[0].split('.')[0]] = (sym[0], os.path.basename(f))
        except Exception:
            pass
loci = []
for r in sorted(rows, key=lambda r: (r['trait'], r['chr'], r['lo'])):
    for L in loci:
        if L['trait'] == r['trait'] and L['chr'] == r['chr'] and (r['lo'] <= L['hi'] + 500000) and (r['hi'] >= L['lo'] - 500000):
            L['hi'] = max(L['hi'], r['hi'])
            L['lo'] = min(L['lo'], r['lo'])
            L['genes'].append(r['gene'].split('.')[0])
            L['pmin'] = min(L['pmin'], r['pmin'])
            L['maxz'] = max(L['maxz'], r['maxz'])
            break
    else:
        loci.append(dict(trait=r['trait'], chr=r['chr'], lo=r['lo'], hi=r['hi'], genes=[r['gene'].split('.')[0]], pmin=r['pmin'], maxz=r['maxz']))
for L in loci:
    f = f"{W}/ref/gwas05/{L['trait']}.chr{L['chr']}.txt"
    mp = 1.0
    mpos = None
    with open(f) as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        ci = {h: i for i, h in enumerate(hdr)}
        for line in fh:
            q = line.split('\t')
            pos = int(q[ci['POS']])
            if L['lo'] - 500000 <= pos <= L['hi'] + 500000:
                try:
                    p = float(q[ci['p.value']])
                except:
                    continue
                if p < mp:
                    mp = p
                    mpos = pos
    L['common_minp_500kb'] = mp
    L['common_minpos'] = mpos
    L['names'] = [names.get(g, ('?', ''))[0] for g in L['genes']]
print('name_table_hits', len(names))
for L in loci:
    print(json.dumps(dict(trait=L['trait'], chr=L['chr'], span_kb=round((L['hi'] - L['lo']) / 1000.0), n_win=len(L['genes']), pmin=L['pmin'], maxz=L['maxz'], common_minp=f"{L['common_minp_500kb']:.1e}", genes=L['genes'][:4], names=L['names'][:4])))
json.dump(loci, open(f'{W}/run_meth/out/gate2_B86_novel_loci.json', 'w'), indent=1)
