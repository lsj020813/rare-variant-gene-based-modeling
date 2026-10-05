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
import json
from pathlib import Path
import sys

sys.path.insert(0, _config_path('${PROJECT_ROOT}/work/run_l3'))
import l3_common as C

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
    args = ap.parse_args()
    C.load_libraries(args.threads)
    np = C.np
    out = Path(args.out)

    maf = {}
    bad = 0
    for line in open(args.maf_tsv):
        k, _, v = line.rstrip('\n').partition('\t')
        try:
            maf[k] = float(v)
        except ValueError:
            bad += 1
    C.require(bool(maf), 'empty MAF table')

    lines, n_slots, n_fallback = [], 0, 0
    phis, betas, prods = [], [], []
    cur_var = None
    for line in open(args.g_phi):
        t = line.split()
        if len(t) < 3:
            continue
        g, kind = t[0], t[1]
        if kind == 'var':
            cur_var = t[2:]
            lines.append(line.rstrip('\n'))
        elif kind == 'anno':
            lines.append(line.rstrip('\n'))
        elif kind == 'weight':
            C.require(cur_var is not None and len(cur_var) == len(t) - 2,
                      f'{g}: weight/var length mismatch')
            ws = []
            for k, wtxt in zip(cur_var, t[2:]):
                p = float(wtxt)
                m = maf.get(k)
                if m is None:
                    n_fallback += 1
                    m = MAF_FALLBACK
                b = beta(m)
                phis.append(p); betas.append(b); prods.append(p * b)
                ws.append(f'{p * b:.8g}')
            n_slots += len(ws)
            lines.append(f'{g} weight ' + ' '.join(ws))
            cur_var = None
        else:
            C.require(False, f'unexpected line kind {kind}')

    tmp = out.with_name(out.name + '.tmp')
    with open(tmp, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
        fh.flush()
    tmp.replace(out)

    a = np.asarray(prods)
    manifest = dict(
        arm='G_phibeta', weight='phi * Beta(MAF; 1, 25)',
        maf_source="bcftools query -f '%CHROM:%POS:%REF:%ALT\\t%INFO/MAF\\n' ref/band_vcf/chr19.band.vcf.gz",
        maf_field='INFO/MAF', beta_formula=BETA_FORMULA, maf_fallback=MAF_FALLBACK,
        provenance='MAF source and Beta formula copied verbatim from run_bwg2/wtest.sh (A2 arm, '
                   'which reproduced SAIGE default Beta weighting to median rel. diff 1.4e-3)',
        maf_table_rows=len(maf), maf_table_unparsed=bad,
        variant_slots=n_slots, maf_fallback_slots=n_fallback,
        phi_source=str(args.g_phi),
        stats=dict(
            phi={k: float(v) for k, v in zip(('mean', 'sd', 'min', 'max'),
                 (np.mean(phis), np.std(phis), np.min(phis), np.max(phis)))},
            beta={k: float(v) for k, v in zip(('mean', 'sd', 'min', 'max'),
                  (np.mean(betas), np.std(betas), np.min(betas), np.max(betas)))},
            product={**{k: float(v) for k, v in zip(('mean', 'sd', 'min', 'max'),
                     (a.mean(), a.std(), a.min(), a.max()))},
                     'q': {f'p{q}': float(np.percentile(a, q)) for q in (1, 25, 50, 75, 99)}}),
        file=dict(path=str(out), bytes=out.stat().st_size, sha256=C.sha_file(out)),
        code_sha256=C.sha_file(__file__))
    C.atomic_json(out.parent / 'gf_phibeta_manifest.json', manifest)
    C.atomic_json(out.parent / 'gf_phibeta_manifest.done',
                  dict(status='complete', slots=n_slots, maf_fallback_slots=n_fallback,
                       sha256=manifest['file']['sha256']))
    C.log(f'PHIBETA_DONE slots={n_slots} maf_fallback={n_fallback} '
          f'w[p1,p50,p99]=({np.percentile(a,1):.5g},{np.percentile(a,50):.5g},{np.percentile(a,99):.5g})')

if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        msg = str(error) if isinstance(error, RuntimeError) and str(error).startswith('GATE FAIL:') \
            else type(error).__name__
        print(f'[mk_gf_phibeta] {msg}', flush=True)
        raise SystemExit(1) from None
