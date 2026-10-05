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


import argparse
import math
from pathlib import Path
import sys

sys.path.insert(0, _config_path('${PROJECT_ROOT}/work/run_l3'))
import l3_common as C

PHI_POWER = 0.5
BETA_FORMULA = 'beta(m) = 25.0 * (1.0 - m) ** 24  # dbeta(m; 1, 25), verbatim from run_bwg2/wtest.sh'
MAF_FALLBACK = 0.005

def beta(m):
    return 25.0 * (1.0 - m) ** 24

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--g-phi', required=True)
    ap.add_argument('--maf-tsv', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--threads', required=True, type=int)
    a = ap.parse_args()
    C.load_libraries(a.threads)
    np = C.np
    out = Path(a.out)

    maf, bad = {}, 0
    for line in open(a.maf_tsv):
        k, _, v = line.rstrip('\n').partition('\t')
        try:
            maf[k] = float(v)
        except ValueError:
            bad += 1
    C.require(bool(maf), 'empty MAF table')

    lines, n_slots, n_fallback = [], 0, 0
    phis, sqrts, betas, prods = [], [], [], []
    cur_var = None
    for line in open(a.g_phi):
        t = line.split()
        if len(t) < 3:
            continue
        g, kind = t[0], t[1]
        if kind in ('var', 'anno'):
            if kind == 'var':
                cur_var = t[2:]
            lines.append(line.rstrip('\n'))
        elif kind == 'weight':
            C.require(cur_var is not None and len(cur_var) == len(t) - 2,
                      f'{g}: weight/var length mismatch')
            ws = []
            for k, wtxt in zip(cur_var, t[2:]):
                p = float(wtxt)
                C.require(p > 0, 'non-positive phi encountered')
                s = p ** PHI_POWER
                m = maf.get(k)
                if m is None:
                    n_fallback += 1
                    m = MAF_FALLBACK
                b = beta(m)
                phis.append(p); sqrts.append(s); betas.append(b); prods.append(s * b)
                ws.append(f'{s * b:.8g}')
            n_slots += len(ws)
            lines.append(f'{g} weight ' + ' '.join(ws))
            cur_var = None
        else:
            C.require(False, f'unexpected line kind {kind}')

    tmp = out.with_name(out.name + '.tmp')
    with open(tmp, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    tmp.replace(out)

    def st(v):
        v = np.asarray(v)
        return dict(mean=float(v.mean()), sd=float(v.std()), min=float(v.min()), max=float(v.max()),
                    q={f'p{q}': float(np.percentile(v, q)) for q in (1, 25, 50, 75, 99)},
                    p99_over_p1=float(np.percentile(v, 99) / np.percentile(v, 1)))
    manifest = dict(
        arm='G_phisqrt', weight='sqrt(phi) * Beta(MAF; 1, 25)', phi_power=PHI_POWER,
        rationale='SKAT uses w^2, so w ~ sqrt(phi) makes the annotation enter linearly in Q; '
                  'pre-registration v108 L3-13. Exponent fixed at 0.5, not tuned.',
        maf_source="bcftools query -f '%CHROM:%POS:%REF:%ALT\\t%INFO/MAF\\n' ref/band_vcf/chr19.band.vcf.gz",
        maf_field='INFO/MAF', beta_formula=BETA_FORMULA, maf_fallback=MAF_FALLBACK,
        provenance='MAF source and Beta formula identical to the phibeta arm (run_bwg2/wtest.sh A2)',
        phi_source=str(a.g_phi), phi_reused_not_recomputed=True,
        maf_table_rows=len(maf), maf_table_unparsed=bad,
        variant_slots=n_slots, maf_fallback_slots=n_fallback,
        stats=dict(phi=st(phis), sqrt_phi=st(sqrts), beta=st(betas), weight=st(prods)),
        file=dict(path=str(out), bytes=out.stat().st_size, sha256=C.sha_file(out)),
        code_sha256=C.sha_file(__file__))
    C.atomic_json(out.parent / 'gf_phisqrt_manifest.json', manifest)
    C.atomic_json(out.parent / 'gf_phisqrt_manifest.done',
                  dict(status='complete', slots=n_slots, maf_fallback_slots=n_fallback,
                       sha256=manifest['file']['sha256']))
    C.log(f'PHISQRT_DONE slots={n_slots} maf_fallback={n_fallback} '
          f'phi_spread_p99/p1={manifest["stats"]["phi"]["p99_over_p1"]:.4g} '
          f'sqrt_spread={manifest["stats"]["sqrt_phi"]["p99_over_p1"]:.4g} '
          f'w[p1,p50,p99]=({np.percentile(prods,1):.5g},{np.percentile(prods,50):.5g},{np.percentile(prods,99):.5g})')

if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        msg = str(error) if isinstance(error, RuntimeError) and str(error).startswith('GATE FAIL:') \
            else type(error).__name__
        print(f'[mk_gf_phisqrt] {msg}', flush=True)
        raise SystemExit(1) from None
