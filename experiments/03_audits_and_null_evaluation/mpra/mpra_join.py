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
import subprocess, gzip, json, sys, collections, os
BCF = _configured('${BCFTOOLS}')
MP = sys.argv[1]
OUT = sys.argv[2]
mp = {}
with gzip.open(MP, 'rt') as f:
    hdr = f.readline().rstrip('\n').split('\t')
    for line in f:
        p = line.rstrip('\n').split('\t')
        d = dict(zip(hdr, p))
        mp[d['chr19'], int(d['pos19'])] = d
res = {}
for band, root, gdir in [('0.1-1%', _configured('${PROJECT_ROOT}/work/ref/band_vcf'), _configured('${PROJECT_ROOT}/work/ref/groupfiles_bwg')), ('1-5%', _configured('${PROJECT_ROOT}/work/ref15/band_vcf'), _configured('${PROJECT_ROOT}/work/ref15/groupfiles_bwg'))]:
    r = collections.Counter()
    win_tested = collections.Counter()
    win_emvar = collections.Counter()
    nwin = 0
    matched = set()
    for N in range(1, 23):
        V = f'{root}/chr{N}.band.vcf.gz'
        if not os.path.exists(V):
            continue
        q = subprocess.Popen([BCF, 'query', '-f', '%CHROM\t%POS\t%REF\t%ALT\t%INFO/MAF\t%INFO/R2\n', V], stdout=subprocess.PIPE, text=True)
        for line in q.stdout:
            c, pos, ref, alt, maf, r2 = line.rstrip('\n').split('\t')
            c = c.replace('chr', '')
            r['band_variants'] += 1
            d = mp.get((c, int(pos)))
            if d is None:
                continue
            r['pos_match'] += 1
            allele_ok = {ref, alt} == {d['ref'], d['alt']}
            if not allele_ok:
                r['pos_match_allele_mismatch'] += 1
                continue
            r['tested'] += 1
            matched.add((c, int(pos), ref, alt))
            if d['emvar_any'] == 'True':
                r['emvar_any'] += 1
            if d['emvar_hepg2'] == 'True':
                r['emvar_hepg2'] += 1
            if d['finemapped_set'] == 'True':
                r['tested_finemapped_set'] += 1
            if float(r2) >= 0.8:
                r['tested_R2ge0.8'] += 1
        q.wait()
        G = f'{gdir}/chr{N}.B_3kb_re2g.txt'
        if os.path.exists(G):
            with open(G) as g:
                for line in g:
                    p = line.split()
                    if len(p) < 3 or p[1] != 'var':
                        continue
                    nwin += 1
                    t = 0
                    e = 0
                    for vid in p[2:]:
                        tok = vid.replace('chr', '').replace('_', ':').split(':')
                        if len(tok) >= 4:
                            key = (tok[0], int(tok[1]), tok[2], tok[3])
                            if key in matched:
                                t += 1
                                d = mp.get((tok[0], int(tok[1])))
                                if d and d['emvar_any'] == 'True':
                                    e += 1
                    win_tested[min(t, 3)] += 1
                    win_emvar[min(e, 3)] += 1
    res[band] = dict(counts=dict(r), windows=nwin, windows_by_n_tested={str(k): v for k, v in sorted(win_tested.items())}, windows_by_n_emvar={str(k): v for k, v in sorted(win_emvar.items())})
    print(band, json.dumps(res[band]), flush=True)
json.dump(res, open(OUT, 'w'), indent=1)
