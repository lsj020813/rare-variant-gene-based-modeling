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
import os
from pathlib import Path
import sys

sys.path.insert(0, _config_path('${PROJECT_ROOT}/work/run_l3'))
import l3_common as C

def parser():
    p = argparse.ArgumentParser()
    p.add_argument('--root', required=True)
    p.add_argument('--cadd-dir', required=True)
    p.add_argument('--cadd-source', required=True, choices=['features', 'extract'])
    p.add_argument('--chrom', required=True)
    p.add_argument('--group-file', required=True)
    p.add_argument('--phi-json', required=True)
    p.add_argument('--v10-script', required=True)
    p.add_argument('--t1-partial-policy', required=True, choices=['blank', 'fail'])
    p.add_argument('--missing-phi-rule', required=True, choices=['intercept-prior', 'fail'],
                   help='intercept-prior = expit(default_intercept); recorded as a star if ever used')
    p.add_argument('--threads', required=True, type=int)
    p.add_argument('--out', required=True)
    return p

def main():
    args = parser().parse_args()
    C.load_libraries(args.threads)
    np = C.np
    from scipy.special import expit
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    v10 = C.import_v10(args.v10_script)
    cache = C.BandCache(args.root, args.cadd_dir, args.cadd_source, [args.chrom], args.t1_partial_policy)
    phi, meta = C.phi_arm(cache, 'learned', v10, args.phi_json, 'band-all')
    C.log(f'cache pairs={len(cache.keys)} genes={cache.ng} phi mean={phi.mean():.6g} '
          f'sd={phi.std():.6g} range=[{phi.min():.6g},{phi.max():.6g}]')

    table = {}
    for i in range(len(cache.keys)):
        table[(cache.genes[i], cache.keys[i])] = float(phi[i])
    C.require(len(table) == len(cache.keys), 'duplicate (gene, variant) pair in the cache')

    with open(args.phi_json) as fh:
        model = json.load(fh)
    fallback = float(expit(model['default_intercept']))

    genes, order = {}, []
    for line in open(args.group_file):
        t = line.split()
        if len(t) < 3:
            continue
        g, kind = t[0], t[1]
        if g not in genes:
            genes[g] = {}
            order.append(g)
        C.require(kind not in genes[g], f'duplicate {kind} line for a gene')
        genes[g][kind] = t[2:]
    C.require(bool(order), 'empty group file')

    n_var = n_miss = 0
    miss_genes = set()
    phis_used = []
    paths = {k: out / f'G_{k}.txt' for k in ('none', 'flat', 'phi')}
    with open(paths['none'], 'w') as f0, open(paths['flat'], 'w') as f1, open(paths['phi'], 'w') as f2:
        for g in order:
            v = genes[g]['var']
            a = genes[g].get('anno', ['all'] * len(v))
            C.require(len(v) == len(a), f'{g}: var/anno length mismatch')
            base = g.split('.')[0]
            ws = []
            for k in v:
                w = table.get((base, k))
                if w is None:
                    n_miss += 1
                    miss_genes.add(base)
                    C.require(args.missing_phi_rule == 'intercept-prior',
                              f'{n_miss} (gene,variant) pairs lack phi; rule=fail')
                    w = fallback
                ws.append(w)
            phis_used.extend(ws)
            n_var += len(v)
            for fh in (f0, f1, f2):
                fh.write(f'{g} var ' + ' '.join(v) + '\n')
                fh.write(f'{g} anno ' + ' '.join(a) + '\n')
            f1.write(f'{g} weight ' + ' '.join(['1.0'] * len(v)) + '\n')
            f2.write(f'{g} weight ' + ' '.join(f'{w:.8g}' for w in ws) + '\n')

    def tokens(path):
        return [tuple(l.split()) for l in open(path) if len(l.split()) >= 3]
    C.require(tokens(paths['none']) == tokens(args.group_file),
              'G_none is not token-identical to the original group file')

    w = np.asarray(phis_used)
    manifest = dict(chrom=args.chrom, group_file=str(args.group_file),
                    genes=len(order), variant_slots=n_var,
                    cache_pairs=int(len(cache.keys)),
                    join_exact=bool(n_var == len(cache.keys) and n_miss == 0),
                    missing_phi_slots=n_miss, missing_phi_genes=len(miss_genes),
                    missing_phi_rule=args.missing_phi_rule, missing_phi_fallback=fallback,
                    phi_scale='probability, unchanged (no logit; D4 unsealed)',
                    phi_stats=dict(mean=float(w.mean()), sd=float(w.std()), min=float(w.min()),
                                   max=float(w.max()),
                                   q={f'p{q}': float(np.percentile(w, q)) for q in (1, 25, 50, 75, 99)}),
                    weight_semantics='REPLACE (verified 2026-09-09 wtest): G_phi drops the MAF Beta weight',
                    phi_meta=meta, t1_partial_blanked_cells=int(cache.t1_partial_blanked),
                    cadd_source=args.cadd_source, cadd_chroms_loaded=len(cache.cadd_available),
                    files={k: dict(path=str(p), bytes=p.stat().st_size, sha256=C.sha_file(p))
                           for k, p in paths.items()},
                    input_fingerprints=cache.fingerprints,
                    code_sha256=C.sha_file(__file__),
                    common_sha256=C.sha_file(_config_path('${PROJECT_ROOT}/work/run_l3/l3_common.py')))
    C.atomic_json(out / 'gf_manifest.json', manifest)
    C.atomic_json(out / 'gf_manifest.done',
                  dict(status='complete', sha256=C.sha_file(out / 'gf_manifest.json'),
                       genes=len(order), variant_slots=n_var, missing_phi_slots=n_miss))
    C.log(f'GF_DONE genes={len(order)} slots={n_var} missing_phi={n_miss} '
          f'phi[p1,p50,p99]=({np.percentile(w,1):.4g},{np.percentile(w,50):.4g},{np.percentile(w,99):.4g})')

if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        msg = str(error) if isinstance(error, RuntimeError) and str(error).startswith('GATE FAIL:') \
            else f'{type(error).__name__}'
        print(f'[mk_gf] {msg}', flush=True)
        raise SystemExit(1) from None
