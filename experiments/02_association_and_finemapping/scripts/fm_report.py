#!/usr/bin/env python3
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import glob, os, gzip, csv, json, sys
W = _config_path('${PROJECT_ROOT}/work/run_ourfm'); G = _config_path('${PROJECT_ROOT}/work/ref/gwas05')
traits = ['tchl', 'htn', 'dm', 'lip']
os.makedirs(f'{W}/report', exist_ok=True)
gw_done = {t: sorted(int(os.path.basename(p).split('.chr')[1].split('.')[0]) for p in glob.glob(f'{G}/{t}.chr*.txt.done')) for t in traits}
chr_all4 = sorted(set.intersection(*[set(v) for v in gw_done.values()]))
status_rows = []; summ = {t: dict(n_chr_regions=0, n_region=0, n_lead=0, n_done=0, n_skip=0, n_cs=0, n_pip01=0, n_pip05=0, n_fm=0,
                                  mass_ge5=0.0, mass_1to5=0.0, mass_0p1to1=0.0, mass_lt0p1=0.0, nvar_ge5=0, nvar_1to5=0, nvar_0p1to1=0, nvar_lt0p1=0,
                                  cs_sizes=[], n_nonconv=0, n_smoke=0) for t in traits}
pip_rows = {t: [] for t in traits}; pip_hdr = None
for t in traits:
    for rf in sorted(glob.glob(f'{W}/regions/{t}.chr*.regions.tsv')):
        chrom = int(rf.split('.chr')[1].split('.')[0])
        with open(rf) as f:
            rows = list(csv.DictReader(f, delimiter='\t'))
        real = [r for r in rows if r['smoke'] == '0']
        summ[t]['n_smoke'] += len(rows) - len(real)
        summ[t]['n_chr_regions'] += 1
        for r in real:
            rid = r['region_id']; d = f'{W}/fm/{t}/{rid}'
            summ[t]['n_region'] += 1; summ[t]['n_lead'] += int(r['n_lead'])
            st = 'pending'
            if os.path.exists(d + '.skip'):
                st = 'skip:' + open(d + '.skip').read().strip(); summ[t]['n_skip'] += 1
            elif os.path.exists(d + '.done') and os.path.exists(d + '.summary.tsv'):
                st = 'done'; summ[t]['n_done'] += 1
                s = list(csv.DictReader(open(d + '.summary.tsv'), delimiter='\t'))[0]
                summ[t]['n_cs'] += int(s['n_cs']); summ[t]['n_pip01'] += int(s['n_pip_ge0.1']); summ[t]['n_pip05'] += int(s['n_pip_ge0.5']); summ[t]['n_fm'] += int(s['n_fm'])
                for k, kk in [('mass_ge5', 'mass_ge5'), ('mass_1to5', 'mass_1to5'), ('mass_0.1to1', 'mass_0p1to1'), ('mass_lt0.1', 'mass_lt0p1')]:
                    summ[t][kk] += float(s[k])
                for k, kk in [('nvar_ge5', 'nvar_ge5'), ('nvar_1to5', 'nvar_1to5'), ('nvar_0.1to1', 'nvar_0p1to1'), ('nvar_lt0.1', 'nvar_lt0p1')]:
                    summ[t][kk] += int(s[k])
                if s['cs_sizes']: summ[t]['cs_sizes'] += [int(x) for x in s['cs_sizes'].split(';')]
                if s['converged'] != 'TRUE': summ[t]['n_nonconv'] += 1
                with open(d + '.pip.tsv') as f:
                    rd = csv.reader(f, delimiter='\t'); h = next(rd)
                    if pip_hdr is None: pip_hdr = h
                    pip_rows[t] += list(rd)
            status_rows.append(dict(trait=t, chr=chrom, region=rid, width_mb=round((int(r['end']) - int(r['start'])) / 1e6, 2), n_lead=r['n_lead'], min_p=r['min_p'], status=st))
    with gzip.open(f'{W}/our_pip/{t}.tsv.gz', 'wt') as f:
        w = csv.writer(f, delimiter='\t'); w.writerow(pip_hdr or ['variant_hg19', 'chr', 'pos', 'a1', 'a2', 'maf', 'z', 'pip', 'cs_id', 'region']); w.writerows(pip_rows[t])
with open(f'{W}/report/regions_status.tsv', 'w') as f:
    w = csv.DictWriter(f, fieldnames=['trait', 'chr', 'region', 'width_mb', 'n_lead', 'min_p', 'status'], delimiter='\t'); w.writeheader(); w.writerows(status_rows)
cols = ['trait', 'gwas_chr_done', 'n_chr_regions', 'n_region', 'n_lead', 'n_done', 'n_skip', 'n_cs', 'cs_size_median', 'cs_size_max', 'n_pip01', 'n_pip05', 'n_fm',
        'nvar_ge5', 'nvar_1to5', 'nvar_0p1to1', 'nvar_lt0p1', 'mass_ge5', 'mass_1to5', 'mass_0p1to1', 'mass_lt0p1', 'n_nonconv', 'n_smoke_excluded']
out = []
tot = dict(n_region=0, n_done=0, n_skip=0, n_cs=0, n_pip01=0)
for t in traits:
    s = summ[t]; cs = sorted(s['cs_sizes'])
    med = cs[len(cs) // 2] if cs else 0
    out.append([t, len(gw_done[t]), s['n_chr_regions'], s['n_region'], s['n_lead'], s['n_done'], s['n_skip'], s['n_cs'], med, max(cs) if cs else 0, s['n_pip01'], s['n_pip05'], s['n_fm'],
                s['nvar_ge5'], s['nvar_1to5'], s['nvar_0p1to1'], s['nvar_lt0p1'], round(s['mass_ge5'], 3), round(s['mass_1to5'], 3), round(s['mass_0p1to1'], 3), round(s['mass_lt0p1'], 3), s['n_nonconv'], s['n_smoke']])
    for k in tot: tot[k] += s[k]
with open(f'{W}/report/summary_by_trait.tsv', 'w') as f:
    w = csv.writer(f, delimiter='\t'); w.writerow(cols); w.writerows(out)
print('\t'.join(cols))
for o in out: print('\t'.join(map(str, o)))
print(json.dumps(dict(chr_all4_gwas_done=chr_all4, n_chr_all4=len(chr_all4), total=tot)))
