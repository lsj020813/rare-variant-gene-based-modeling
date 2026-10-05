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


import sys
trait, chrom = sys.argv[1], sys.argv[2]
smoke_top = '--smoke-top' in sys.argv[3:]
gw = _config_path(f'${{PROJECT_ROOT}}/work/ref/gwas05/{trait}.chr{chrom}.txt')
out = _config_path(f'${{PROJECT_ROOT}}/work/run_ourfm/regions/{trait}.chr{chrom}.regions.tsv')
P_THR = 5e-8; WIN = 1500000; MHC = (25000000, 36000000)
leads = []; n = 0; minp = 1.0; minpos = None; mhc_lead_drop = 0
with open(gw) as f:
    hdr = f.readline().rstrip('\n').split('\t')
    ip = hdr.index('p.value'); ipos = hdr.index('POS')
    for line in f:
        s = line.rstrip('\n').split('\t'); n += 1
        try:
            p = float(s[ip]); pos = int(s[ipos])
        except ValueError:
            continue
        if p < minp: minp = p; minpos = pos
        if p < P_THR:
            if chrom == '6' and MHC[0] <= pos <= MHC[1]:
                mhc_lead_drop += 1; continue
            leads.append((pos, p))
smoke = False
if not leads and smoke_top and minpos is not None:
    leads = [(minpos, minp)]; smoke = True
leads.sort()
regions = []
for pos, p in leads:
    s, e = max(1, pos - WIN), pos + WIN
    if regions and s <= regions[-1][1]:
        regions[-1][1] = max(regions[-1][1], e); regions[-1][2] += 1; regions[-1][3] = min(regions[-1][3], p)
    else:
        regions.append([s, e, 1, p])
mhc_overlap = sum(1 for r in regions if chrom == '6' and r[0] <= MHC[1] and r[1] >= MHC[0])
with open(out, 'w') as o:
    o.write('region_id\tchr\tstart\tend\tn_lead\tmin_p\tsmoke\n')
    for i, r in enumerate(regions, 1):
        o.write(f'{trait}_chr{chrom}_r{i}\t{chrom}\t{r[0]}\t{r[1]}\t{r[2]}\t{r[3]:.3g}\t{int(smoke)}\n')
widths = [ (r[1]-r[0])/1e6 for r in regions ]
print(f'[regions {trait} chr{chrom}] n_var={n} n_lead={len(leads)} n_region={len(regions)} '
      f'mhc_lead_drop={mhc_lead_drop} mhc_region_overlap={mhc_overlap} minp={minp:.3g} '
      f'width_Mb_max={max(widths) if widths else 0:.1f} smoke={int(smoke)}')
