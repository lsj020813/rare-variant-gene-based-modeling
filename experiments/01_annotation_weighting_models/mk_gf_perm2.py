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
from pathlib import Path
import sys

sys.path.insert(0, _config_path('${PROJECT_ROOT}/work/run_l3'))
import l3_common as C

BETA_FORMULA = 'beta(m) = 25.0 * (1.0 - m) ** 24  # dbeta(m; 1, 25), verbatim from run_bwg2/wtest.sh'
MAF_FALLBACK = 0.005
SQRT_POWER = 0.5

def beta(m):
    return 25.0 * (1.0 - m) ** 24

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--g-phi', required=True)
    ap.add_argument('--maf-tsv', required=True)
    ap.add_argument('--base-form', required=True, choices=['phi', 'sqrt_phi'])
    ap.add_argument('--shuffle-unit', required=True, choices=['within_gene', 'global'])
    ap.add_argument('--seed', required=True, type=int)
    ap.add_argument('--out', required=True)
    ap.add_argument('--threads', required=True, type=int)
    a = ap.parse_args()
    C.load_libraries(a.threads)
    np = C.np
    out = Path(a.out)
    rng = np.random.default_rng(a.seed)

    maf, bad = {}, 0
    for line in open(a.maf_tsv):
        k, _, v = line.rstrip('\n').partition('\t')
        try:
            maf[k] = float(v)
        except ValueError:
            bad += 1
    C.require(bool(maf), 'empty MAF table')

    order, block = [], {}
    cur = None
    for line in open(a.g_phi):
        t = line.split()
        if len(t) < 3:
            continue
        g, kind = t[0], t[1]
        if g not in block:
            block[g] = {}
            order.append(g)
        if kind == 'var':
            block[g]['var'] = t[2:]
            cur = g
        elif kind == 'anno':
            block[g]['anno'] = t[2:]
        elif kind == 'weight':
            block[g]['phi'] = np.array([float(x) for x in t[2:]], dtype=float)
        else:
            C.require(False, f'unexpected line kind {kind}')
    C.require(all({'var', 'anno', 'phi'} <= set(block[g]) for g in order), 'incomplete gene block')
    for g in order:
        C.require(len(block[g]['var']) == len(block[g]['phi']) == len(block[g]['anno']),
                  f'{g}: length mismatch')

    base_of = (lambda p: p) if a.base_form == 'phi' else (lambda p: p ** SQRT_POWER)

    moved = same = 0
    if a.shuffle_unit == 'within_gene':
        for g in order:
            b = base_of(block[g]['phi'])
            perm = rng.permutation(len(b))
            moved += int((perm != np.arange(len(b))).sum())
            same += int((perm == np.arange(len(b))).sum())
            block[g]['base'] = b[perm]
    else:
        flat = np.concatenate([base_of(block[g]['phi']) for g in order])
        perm = rng.permutation(len(flat))
        moved = int((perm != np.arange(len(flat))).sum())
        same = int((perm == np.arange(len(flat))).sum())
        shuffled = flat[perm]
        i = 0
        for g in order:
            n = len(block[g]['phi'])
            block[g]['base'] = shuffled[i:i + n]
            i += n
        C.require(i == len(flat), 'global reassignment length mismatch')

    lines, n_slots, n_fallback = [], 0, 0
    bases, prods = [], []
    for g in order:
        v, an, b = block[g]['var'], block[g]['anno'], block[g]['base']
        ws = []
        for k, x in zip(v, b):
            m = maf.get(k)
            if m is None:
                n_fallback += 1
                m = MAF_FALLBACK
            w = float(x) * beta(m)
            bases.append(float(x)); prods.append(w)
            ws.append(f'{w:.8g}')
        n_slots += len(ws)
        lines.append(f'{g} var ' + ' '.join(v))
        lines.append(f'{g} anno ' + ' '.join(an))
        lines.append(f'{g} weight ' + ' '.join(ws))
    tmp = out.with_name(out.name + '.tmp')
    with open(tmp, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    tmp.replace(out)

    if a.shuffle_unit == 'within_gene':
        for g in order:
            C.require(np.allclose(np.sort(block[g]['base']), np.sort(base_of(block[g]['phi'])),
                                  rtol=0, atol=0), f'{g}: within-gene multiset not preserved')

    def st(v):
        v = np.asarray(v)
        return dict(mean=float(v.mean()), sd=float(v.std()), min=float(v.min()), max=float(v.max()),
                    q={f'p{q}': float(np.percentile(v, q)) for q in (1, 25, 50, 75, 99)})
    manifest = dict(
        arm=out.stem,
        weight=('permuted(phi) * Beta(MAF; 1, 25)' if a.base_form == 'phi'
                else 'permuted(sqrt(phi)) * Beta(MAF; 1, 25)'),
        base_form=a.base_form, sqrt_power=(SQRT_POWER if a.base_form == 'sqrt_phi' else None),
        shuffle_unit=a.shuffle_unit, seed=a.seed,
        reference_arm=('phibeta' if a.base_form == 'phi' else 'phisqrt'),
        design_note=('base value distribution preserved exactly; only the base-to-variant '
                     'assignment is randomised. Beta(MAF) stays attached to its own variant.'),
        within_gene_multiset_preserved=(a.shuffle_unit == 'within_gene'),
        maf_source="bcftools query -f '%CHROM:%POS:%REF:%ALT\\t%INFO/MAF\\n' ref/band_vcf/chr19.band.vcf.gz",
        maf_field='INFO/MAF', beta_formula=BETA_FORMULA, maf_fallback=MAF_FALLBACK,
        phi_source=str(a.g_phi), phi_source_sha256=C.sha_file(a.g_phi),
        phi_reused_not_recomputed=True,
        genes=len(order), variant_slots=n_slots, maf_fallback_slots=n_fallback,
        slots_moved_by_permutation=moved, slots_left_in_place=same,
        maf_table_rows=len(maf), maf_table_unparsed=bad,
        stats=dict(permuted_base=st(bases), weight=st(prods)),
        file=dict(path=str(out), bytes=out.stat().st_size, sha256=C.sha_file(out)),
        code_sha256=C.sha_file(__file__))
    C.atomic_json(out.parent / f'gf_{out.stem}_manifest.json', manifest)
    C.atomic_json(out.parent / f'gf_{out.stem}_manifest.done',
                  dict(status='complete', arm=out.stem, base_form=a.base_form,
                       shuffle_unit=a.shuffle_unit, seed=a.seed, variant_slots=n_slots,
                       sha256=manifest['file']['sha256']))
    C.log(f'{out.stem}_DONE base={a.base_form} unit={a.shuffle_unit} seed={a.seed} '
          f'genes={len(order)} slots={n_slots} moved={moved} in_place={same} '
          f'maf_fallback={n_fallback}')

if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        msg = str(error) if isinstance(error, RuntimeError) and str(error).startswith('GATE FAIL:') \
            else type(error).__name__
        print(f'[mk_gf_perm2] {msg}', flush=True)
        raise SystemExit(1) from None
