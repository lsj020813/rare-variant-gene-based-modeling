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


import argparse, sys
from pathlib import Path
sys.path.insert(0, _config_path('${PROJECT_ROOT}/work/run_l3'))
import l3_common as C

ALLOWED = (0.1, 0.25, 0.5)
MAF_FALLBACK = 0.005
def beta(m): return 25.0 * (1.0 - m) ** 24

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--alpha', required=True, type=float)
    ap.add_argument('--g-phi', required=True); ap.add_argument('--maf-tsv', required=True)
    ap.add_argument('--out', required=True); ap.add_argument('--threads', required=True, type=int)
    a = ap.parse_args()
    C.require(any(abs(a.alpha - x) < 1e-12 for x in ALLOWED), f'GATE FAIL: alpha {a.alpha} not in pre-registered grid {ALLOWED}')
    C.load_libraries(a.threads); np = C.np
    out = Path(a.out); tag = out.stem.replace('G_', '')
    maf, bad = {}, 0
    for line in open(a.maf_tsv):
        k, _, v = line.rstrip('\n').partition('\t')
        try: maf[k] = float(v)
        except ValueError: bad += 1
    C.require(bool(maf), 'empty MAF table')
    lines, n_slots, n_fb, cur_var = [], 0, 0, None
    rs, ws_all, neg = [], [], 0
    for line in open(a.g_phi):
        t = line.split()
        if len(t) < 3: continue
        g, kind = t[0], t[1]
        if kind in ('var', 'anno'):
            if kind == 'var': cur_var = t[2:]
            lines.append(line.rstrip('\n'))
        elif kind == 'weight':
            C.require(cur_var is not None and len(cur_var) == len(t) - 2, f'{g}: weight/var length mismatch')
            phi = np.array([float(x) for x in t[2:]], dtype=float)
            C.require(bool((phi > 0).all()), 'non-positive phi encountered')
            r = phi / phi.mean()
            ws = []
            for k, rv in zip(cur_var, r):
                m = maf.get(k)
                if m is None: n_fb += 1; m = MAF_FALLBACK
                w = beta(m) * (1.0 + a.alpha * (rv - 1.0))
                if w <= 0: neg += 1
                rs.append(rv); ws_all.append(w); ws.append(f'{w:.8g}')
            n_slots += len(ws); lines.append(f'{g} weight ' + ' '.join(ws)); cur_var = None
        else:
            C.require(False, f'unexpected line kind {kind}')
    C.require(neg == 0, f'GATE FAIL: {neg} non-positive weights (alpha too large for min r_v)')
    tmp = out.with_name(out.name + '.tmp')
    with open(tmp, 'w') as fh: fh.write('\n'.join(lines) + '\n')
    tmp.replace(out)
    def st(v):
        v = np.asarray(v)
        return dict(mean=float(v.mean()), sd=float(v.std()), min=float(v.min()), max=float(v.max()),
                    p99_over_p1=float(np.percentile(v, 99) / np.percentile(v, 1)))
    man = dict(arm=out.stem, alpha=a.alpha, weight='Beta(MAF;1,25) * [1 + alpha*(r_v - 1)], r_v = phi_v / mean_gene(phi)',
               grid_allowed=list(ALLOWED), preregistration='L3-27 (2026-09-10, before results)',
               phi_source=str(a.g_phi), phi_reused_not_recomputed=True, maf_field='INFO/MAF',
               beta_formula='25.0*(1-m)**24', maf_fallback=MAF_FALLBACK, maf_table_rows=len(maf), maf_table_unparsed=bad,
               variant_slots=n_slots, maf_fallback_slots=n_fb, nonpositive_weights=neg,
               stats=dict(r=st(rs), weight=st(ws_all)),
               file=dict(path=str(out), bytes=out.stat().st_size, sha256=C.sha_file(out)), code_sha256=C.sha_file(__file__))
    C.atomic_json(out.parent / f'gf_{tag}_manifest.json', man)
    C.atomic_json(out.parent / f'gf_{tag}_manifest.done', dict(status='complete', slots=n_slots, sha256=man['file']['sha256']))
    C.log(f'ALPHA_DONE {tag} alpha={a.alpha} slots={n_slots} maf_fallback={n_fb} r_min={min(rs):.4g} r_max={max(rs):.4g} '
          f'w_p99/p1={man["stats"]["weight"]["p99_over_p1"]:.4g}')

if __name__ == '__main__':
    try: main()
    except Exception as e:
        print(f'[mk_gf_alpha] {str(e) if str(e).startswith("GATE FAIL:") else type(e).__name__}', flush=True)
        raise SystemExit(1) from None
