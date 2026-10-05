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
import subprocess, json, sys, collections
B = _configured('${BCFTOOLS}')
out = {}
for N in range(1, 23):
    G = _configured(f'${{PROJECT_ROOT}}/work/ref15/groupfiles_bwg/chr{N}.B_3kb_re2g.txt')
    V = _configured(f'${{PROJECT_ROOT}}/work/ref15/band_vcf/chr{N}.band.vcf.gz')
    typed = set()
    r2 = {}
    p = subprocess.Popen([B, 'query', '-f', '%CHROM:%POS:%REF:%ALT\t%INFO/TYPED\t%INFO/R2\n', V], stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        vid, t, r = line.rstrip('\n').split('\t')
        if t != '.':
            typed.add(vid)
        r2[vid] = float(r) if r != '.' else float('nan')
    p.wait()
    hist = collections.Counter()
    hist9 = collections.Counter()
    nwin = 0
    M_tot = 0
    for line in open(G):
        q = line.split()
        if len(q) > 2 and q[1] == 'var':
            nwin += 1
            vids = q[2:]
            M_tot += len(vids)
            nt = sum((1 for v in vids if v in typed))
            hist[min(nt, 5)] += 1
            n9 = sum((1 for v in vids if v in typed or r2.get(v, 0) >= 0.9))
            hist9[min(n9, 10)] += 1
    out[N] = dict(windows=nwin, variants_in_windows=M_tot, typed_total=len(typed), win_by_ntyped={str(k): hist[k] for k in sorted(hist)}, win_by_n_typed_or_R2ge0_9={str(k): hist9[k] for k in sorted(hist9)})
    print(f'chr{N}\twindows={nwin}\ttyped={len(typed)}\twin_typed>=2={sum((v for k, v in hist.items() if k >= 2))}\twin_(typed|R2>=.9)>=2={sum((v for k, v in hist9.items() if k >= 2))}', flush=True)
json.dump(out, open(sys.argv[1], 'w'), indent=1)
print('DONE')
