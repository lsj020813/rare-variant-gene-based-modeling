import os as _cfg_os
import re as _cfg_re
import math as _cfg_math

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)

def _config_number(name, cast=int, positive=False):
    value = cast(_cfg_os.environ[name])
    if not _cfg_math.isfinite(value) or (positive and value <= 0):
        raise ValueError(name + " must be finite" + (" and positive" if positive else ""))
    return value

import os as _os
N_SAMPLES = _config_number("N_SAMPLES", int, True)
import os
import sys
import json
import numpy as np
from cyvcf2 import VCF

import f3f6_common as C

GENO = C.OUT + '/geno'
NSUB = 5000
NTOT = N_SAMPLES
KTOP = 21

def sub_idx():
    return np.sort(np.random.default_rng(C.SEED).choice(NTOT, NSUB,
                                                        replace=False))

def main():
    ch = sys.argv[1]
    os.makedirs(GENO, exist_ok=True)
    d = C.load()
    sel = np.load(C.OUT + '/f46_domains.npy')
    lab = d['mod_spectral_eucl']
    di = d['dom_idx']
    cols = sub_idx()
    vcf = VCF(_config_path('${PROJECT_ROOT}/work/ref/orig_index/chr%s.vcf.gz') % ch)
    assert len(vcf.samples) == NTOT, 'sample 수 불일치 %d' % len(vcf.samples)
    dom_here = [g for g in sel
                if d['dom_chrom'][g].replace('chr', '') == str(ch)]
    rows_f6 = []
    rep = []
    for g in dom_here:
        rid = np.where(di == g)[0]
        p19 = d['pos19'][rid]
        want = {}
        for j, r in enumerate(rid):
            want[d['key'][r]] = j
        M = len(rid)
        D = np.full((NSUB, M), np.nan, dtype=np.float32)
        got = np.zeros(M, dtype=bool)
        sp = np.sort(np.unique(p19))
        iv = []
        a0 = int(sp[0])
        b0 = int(sp[0])
        for x in sp[1:]:
            x = int(x)
            if x - b0 <= 100:
                b0 = x
            else:
                iv.append((a0, b0))
                a0 = x
                b0 = x
        iv.append((a0, b0))
        recs = 0
        for (qa, qb) in iv:
          for v in vcf('%s:%d-%d' % (ch, qa, qb)):
            recs += 1
            for a in v.ALT:
                k = 'chr%s:%d:%s:%s' % (ch, v.POS, v.REF, a)
                j = want.get(k)
                if j is None or got[j]:
                    continue
                try:
                    ds = np.asarray(v.format('DS'), dtype=np.float32).ravel()
                except Exception:
                    ds = None
                if ds is None or ds.shape[0] != NTOT:
                    continue
                x = ds[cols]
                bad = ~np.isfinite(x) | (x < 0)
                if bad.all():
                    continue
                if bad.any():
                    x = np.where(bad, np.nanmean(x[~bad]), x)
                D[:, j] = x
                got[j] = True
        mrate = float(got.mean())
        if got.sum() < 2:
            rep.append(dict(domain_idx=int(g), n_var=M, match_rate=mrate,
                            status='too_few_matched'))
            continue
        keep = np.where(got)[0]
        Dk = D[:, keep]
        sd = Dk.std(axis=0)
        ok = sd > 1e-8
        Z = np.where(ok, (Dk - Dk.mean(0)) / np.where(ok, sd, 1.0), np.nan)
        Zf = np.nan_to_num(Z, nan=0.0)
        R = (Zf.T @ Zf) / float(NSUB)
        R2 = R ** 2
        R2[~ok, :] = np.nan
        R2[:, ~ok] = np.nan
        np.savez_compressed(
            GENO + '/dom%06d.npz' % g,
            local_rows=keep.astype(np.int32), r2=R2.astype(np.float32),
            corr=R.astype(np.float32), var_ds=Dk.var(axis=0),
            ac=Dk.sum(axis=0), carrier_n=(Dk >= 0.5).sum(axis=0),
            nonzero_var=ok, match_rate=np.array([mrate]),
            n_sub=np.array([NSUB]))
        Ck = (Dk >= 0.5)
        lk = lab[rid[keep]]
        for u in np.unique(lk):
            if u < 0:
                continue
            mm = lk == u
            nv = int(mm.sum())
            Cm = Ck[:, mm]
            kper = Cm.sum(axis=1)
            burden = Dk[:, mm].sum(axis=1)
            car = int((kper >= 1).sum())
            hist = np.bincount(np.minimum(kper, KTOP - 1),
                               minlength=KTOP).tolist()
            if nv >= 2:
                cc = Cm.astype(np.float32)
                inter = cc.T @ cc
                cn = Cm.sum(axis=0).astype(np.float64)
                un = cn[:, None] + cn[None, :] - inter
                iu = np.triu_indices(nv, 1)
                with np.errstate(invalid='ignore', divide='ignore'):
                    jac = np.where(un[iu] > 0, inter[iu] / un[iu], np.nan)
                cojac = float(np.nanmean(jac)) if np.isfinite(jac).any() \
                    else float('nan')
                pair_both = float(np.nanmean(inter[iu]))
            else:
                cojac = float('nan')
                pair_both = float('nan')
            rows_f6.append([int(g), int(d['dom_forced'][g]), int(u), nv,
                            car, car / float(NSUB),
                            car * NTOT / float(NSUB),
                            float(burden.var()),
                            float(Dk[:, mm].var(axis=0).sum()),
                            float(kper[kper >= 1].mean()) if car else 0.0,
                            int(kper.max()), cojac, pair_both,
                            int((kper >= 2).sum()), int((kper >= 3).sum()),
                            mrate] + hist)
        rep.append(dict(domain_idx=int(g), n_var=M, n_matched=int(got.sum()),
                        match_rate=mrate, n_records_parsed=recs,
                        n_intervals=len(iv), status='ok'))
    hdr = ['domain_idx', 'forced', 'module', 'n_var', 'n_carriers',
           'carrier_frac_sub', 'carriers_extrap_full_cohort', 'burden_var',
           'sum_var_ds', 'mean_K_per_carrier', 'max_K', 'cocarry_jaccard',
           'mean_pair_both_carriers', 'n_ge2_var', 'n_ge3_var',
           'variant_match_rate'] + ['K%d' % i for i in range(KTOP)]
    with open(C.OUT + '/f6_raw_chr%s.csv' % ch, 'w') as fh:
        fh.write(','.join(hdr) + '\n')
        for r in rows_f6:
            fh.write(','.join([str(x) for x in r]) + '\n')
    with open(C.OUT + '/f46_geno_chr%s.json' % ch, 'w') as fh:
        json.dump(rep, fh)
    with open(C.OUT + '/f46_geno_chr%s.done' % ch, 'w') as fh:
        fh.write('%d domains\n' % len(dom_here))

if __name__ == '__main__':
    main()
