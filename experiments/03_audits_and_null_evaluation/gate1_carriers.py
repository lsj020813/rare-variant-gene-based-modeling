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
import os, sys, subprocess, json, collections, time
BCF = _configured('${BCFTOOLS}')
ROOT = _configured('${COHORT_DATA}')
V = sys.argv[1]
OUT = sys.argv[2]
t0 = time.time()
vcf_ids = subprocess.run([BCF, 'query', '-l', V], capture_output=True, text=True, check=True).stdout.split()
stems = set((f.rsplit('_', 1)[0] for f in os.listdir(_required('HEXA_IDAT_DIR')) if f.endswith('_Grn.idat')))
sub = [s for s in vcf_ids if s in stems]
n = len(sub)
assert n > 0
rfd, wfd = os.pipe()
os.write(wfd, ('\n'.join(sub) + '\n').encode())
os.close(wfd)
view = subprocess.Popen([BCF, 'view', '-I', '-S', '/dev/fd/%d' % rfd, '-Ou', V], stdout=subprocess.PIPE, pass_fds=(rfd,))
q = subprocess.Popen([BCF, 'query', '-f', '%INFO/AC\t%INFO/AN[\t%GT]\n'], stdin=view.stdout, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
view.stdout.close()
bins = [(0, 0, '0'), (1, 4, '1-4'), (5, 9, '5-9'), (10, 19, '10-19'), (20, 10 ** 9, '20+')]

def b(c):
    for lo, hi, l in bins:
        if lo <= c <= hi:
            return l
tot = 0
band = 0
hist = collections.Counter()
miss_tot = 0
maf_bins = collections.Counter()
car_sum = 0
carriers_by_mafbin = collections.defaultdict(list)
for line in q.stdout:
    p = line.rstrip('\n').split('\t')
    tot += 1
    try:
        ac = float(p[0])
        an = float(p[1])
    except ValueError:
        continue
    if an <= 0:
        continue
    af = ac / an
    maf = min(af, 1 - af)
    mb = '<1%' if maf < 0.01 else '1-2%' if maf < 0.02 else '2-5%' if maf <= 0.05 else '>5%'
    maf_bins[mb] += 1
    if not 0.01 <= maf <= 0.05:
        continue
    band += 1
    gts = p[2:]
    car = 0
    miss = 0
    for g in gts:
        if g[0] == '.':
            miss += 1
        elif g != '0/0' and g != '0|0':
            car += 1
    if af > 0.5:
        car = sum((1 for g in gts if g[0] != '.' and '0' in g))
    hist[b(car)] += 1
    miss_tot += miss
    car_sum += car
    carriers_by_mafbin[mb].append(car)
q.wait()
view.wait()

def summ(v):
    v = sorted(v)
    return {'n': len(v), 'median': v[len(v) // 2] if v else None, 'p10': v[len(v) // 10] if v else None, 'p90': v[9 * len(v) // 10] if v else None}
out = {'n_subset': n, 'vcf_total_records': tot, 'maf_bins_cohort': dict(maf_bins), 'band_1_5pct_variants': band, 'carrier_hist': dict(hist), 'mean_carriers_in_band': round(car_sum / band, 2) if band else None, 'missing_rate_in_band': round(miss_tot / (band * n), 5) if band else None, 'carriers_by_mafbin': {k: summ(v) for k, v in carriers_by_mafbin.items()}, 'bcftools_exit': [view.returncode, q.returncode], 'elapsed_s': round(time.time() - t0, 1)}
json.dump(out, open(OUT, 'w'), indent=1)
print(json.dumps(out, indent=1))
