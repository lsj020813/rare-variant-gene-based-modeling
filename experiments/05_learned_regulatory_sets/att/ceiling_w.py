import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import os
import sys
import json
import numpy as np

sys.path.insert(0, _config_path('${PROJECT_ROOT}/work/fset/att'))
import att_common as C

TRAIT = os.environ.get('TRAIT', 'tchl')
C.PHENO = _config_path('${PROJECT_ROOT}/work/ref/pheno_v3/%s_v3.tsv') % TRAIT
import att_fit_lamfix as F

WGRID = [float(x) for x in
         os.environ.get('WGRID', '0,0.25,0.5,0.75,1').split(',')]

def main():
    Z = np.load(C.ATT + '/att_labels.npz', allow_pickle=True)
    gok = Z['dom_ok']
    fold_all = F.person_folds()
    idx, y, Xc = F.load_pheno()
    fold = fold_all[idx]
    print('[%s] N=%d' % (TRAIT, len(y)), flush=True)

    totlab = F.load_col(C.C_TOT_LAB, idx)
    maj = F.load_col(C.C_MAJ0, idx)
    mino = totlab - maj
    tot = F.load_col(C.C_TOT_ALL, idx)

    def rowcorr(A, B):
        a = A - A.mean(1, keepdims=True)
        b = B - B.mean(1, keepdims=True)
        na = np.sqrt((a * a).sum(1))
        nb = np.sqrt((b * b).sum(1))
        ok = (na > 1e-9) & (nb > 1e-9)
        r = np.full(len(A), np.nan)
        r[ok] = ((a[ok] * b[ok]).sum(1) / (na[ok] * nb[ok]))
        return r

    r_mm = rowcorr(maj, mino)

    Zs = {}
    for w in WGRID:
        S = F.arm_fixed_w(maj, mino, np.full(C.NDOM, w, dtype=np.float64),
                          y, fold)
        _, z = F.zstats(S, y, Xc)
        Zs[w] = np.abs(z)
        print('  w=%.2f  중앙 |z| %.3f' % (w, np.nanmedian(np.abs(z))),
              flush=True)

    M = np.vstack([Zs[w] for w in WGRID])
    zmax = np.nanmax(M, 0)
    zmin = np.nanmin(M, 0)
    _, zflat = F.zstats(F.arm_single(tot, y, fold), y, Xc)
    zflat = np.abs(zflat)
    span = zmax - zmin
    gain = zmax - zflat

    m = gok & np.isfinite(span) & np.isfinite(zflat)
    top = zflat >= np.nanpercentile(zflat[m], 95)

    def q(v, mk):
        v = v[mk & np.isfinite(v)]
        return dict(median=float(np.median(v)), p90=float(np.percentile(v, 90)),
                    p99=float(np.percentile(v, 99)), max=float(v.max()))

    res = dict(
        trait=TRAIT, n_individuals=int(len(y)), n_domains=int(m.sum()),
        w_grid=WGRID,
        corr_major_minor=q(r_mm, m),
        zspan_all=q(span, m),
        zspan_top5pct=q(span, m & top),
        gain_over_flat_all=q(gain, m),
        gain_over_flat_top5pct=q(gain, m & top),
        frac_span_gt_0p5=float(np.nanmean(span[m] > 0.5)),
        frac_gain_gt_0p5=float(np.nanmean(gain[m] > 0.5)),
        note=('span = w 격자에서 |z| 의 최대-최소. 어떤 가중 방법이든 '
              '얻을 수 있는 이득의 상한(oracle w 를 안다고 가정).'))
    print()
    print(json.dumps(res, indent=1, ensure_ascii=False))
    out = C.ATT + '/ceiling_w_%s.json' % TRAIT
    with open(out, 'w') as fh:
        json.dump(res, fh, indent=1, ensure_ascii=False)
    print('->', out)

if __name__ == '__main__':
    main()
