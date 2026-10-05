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
import hashlib
from pathlib import Path
import sys

sys.path.insert(0, _config_path('${PROJECT_ROOT}/work/run_l3'))
import l3_common as C

BETA_FORMULA = 'beta(m) = 25.0 * (1.0 - m) ** 24  # dbeta(m; 1, 25), verbatim from run_bwg2/wtest.sh'
MAF_FALLBACK = 0.005
SHUFFLE_UNIT = 'within-gene'

def beta(m):
    return 25.0 * (1.0 - m) ** 24

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--g-phi', required=True)
    ap.add_argument('--maf-tsv', required=True)
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

    lines, n_slots, n_fallback, n_genes = [], 0, 0, 0
    phis, prods = [], []
    moved = same = 0
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
            p = np.array([float(x) for x in t[2:]], dtype=float)
            perm = rng.permutation(len(p))
            pp = p[perm]
            moved += int((perm != np.arange(len(p))).sum())
            same += int((perm == np.arange(len(p))).sum())
            ws = []
            for k, w in zip(cur_var, pp):
                m = maf.get(k)
                if m is None:
                    n_fallback += 1
                    m = MAF_FALLBACK
                v = w * beta(m)
                phis.append(float(w)); prods.append(v)
                ws.append(f'{v:.8g}')
            n_slots += len(ws); n_genes += 1
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
                    q={f'p{q}': float(np.percentile(v, q)) for q in (1, 25, 50, 75, 99)})
    manifest = dict(
        arm=out.stem, weight='permuted(phi) * Beta(MAF; 1, 25)', seed=a.seed,
        shuffle_unit=SHUFFLE_UNIT,
        design_note='phi value distribution preserved exactly; only the phi-to-variant assignment '
                    'is randomised, within gene. Between-gene mean phi is preserved.',
        maf_source="bcftools query -f '%CHROM:%POS:%REF:%ALT\\t%INFO/MAF\\n' ref/band_vcf/chr19.band.vcf.gz",
        maf_field='INFO/MAF', beta_formula=BETA_FORMULA, maf_fallback=MAF_FALLBACK,
        phi_source=str(a.g_phi), phi_source_sha256=C.sha_file(a.g_phi),
        phi_reused_not_recomputed=True,
        genes=n_genes, variant_slots=n_slots, maf_fallback_slots=n_fallback,
        slots_moved_by_permutation=moved, slots_left_in_place=same,
        maf_table_rows=len(maf), maf_table_unparsed=bad,
        stats=dict(permuted_phi=st(phis), weight=st(prods)),
        file=dict(path=str(out), bytes=out.stat().st_size, sha256=C.sha_file(out)),
        code_sha256=C.sha_file(__file__))
    C.atomic_json(out.parent / f'gf_{out.stem}_manifest.json', manifest)
    C.atomic_json(out.parent / f'gf_{out.stem}_manifest.done',
                  dict(status='complete', seed=a.seed, slots=n_slots,
                       sha256=manifest['file']['sha256']))
    C.log(f'{out.stem.upper()}_DONE seed={a.seed} genes={n_genes} slots={n_slots} '
          f'moved={moved} in_place={same} maf_fallback={n_fallback}')

if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        msg = str(error) if isinstance(error, RuntimeError) and str(error).startswith('GATE FAIL:') \
            else type(error).__name__
        print(f'[mk_gf_perm] {msg}', flush=True)
        raise SystemExit(1) from None
