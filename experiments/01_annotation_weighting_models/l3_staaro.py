#!/usr/bin/env python3
import argparse
import json
import math
import os
from pathlib import Path
import resource

import l3_common as C

def parser():
    p = argparse.ArgumentParser(description='layer-3 STAAR-O control (C3)')
    p.add_argument('--root', required=True)
    p.add_argument('--cadd-dir', required=True)
    p.add_argument('--cadd-source', required=True, choices=['features', 'extract'])
    p.add_argument('--chrom', required=True, help="comma list, or 'all'")
    p.add_argument('--traits', required=True)
    p.add_argument('--resid-dir', required=True)
    p.add_argument('--ecdf-reference', required=True, choices=['band-all', 'band-chrom', 'band-cadd-available'])
    p.add_argument('--t1-partial-policy', required=True, choices=['blank', 'fail'])
    p.add_argument('--beta-weights', required=True,
                   help="semicolon list of a,b pairs, e.g. '1,1;1,25' (the paper's two weights)")
    p.add_argument('--acat-set', required=True,
                   choices=['annot-only-burden', 'annot-x-weight', 'annot-x-weight-x-test'],
                   help='which p-values enter the equal-weight ACAT; the paper combines all of them')
    p.add_argument('--staar-null', required=True, choices=['residual-analytic', 'permutation'])
    p.add_argument('--skat-tail', required=True, choices=['imhof', 'liu'],
                   help='mixture-of-chi-square tail; imhof falls back to liu on integration failure')
    p.add_argument('--perm-b', required=True, type=int)
    p.add_argument('--perm-seed', required=True, type=int)
    p.add_argument('--max-variants-per-gene', required=True, type=int)
    p.add_argument('--min-n', required=True, type=int)
    p.add_argument('--threads', required=True, type=int)
    p.add_argument('--memory-gb', required=True, type=float)
    p.add_argument('--min-avail-gb', required=True, type=float)
    p.add_argument('--resid-shuffle-seed', required=True)
    p.add_argument('--perm-strata', required=True)
    p.add_argument('--allow-real-residuals', action='store_true')
    p.add_argument('--out', required=True)
    p.add_argument('--tag', required=True)
    return p

def guard(args):
    with open('/proc/meminfo') as fh:
        info = {k.strip(): v for k, v in (l.split(':') for l in fh)}
    avail = int(info['MemAvailable'].split()[0]) / (1024 ** 2)
    C.require(avail >= args.min_avail_gb, f'MemAvailable {avail:.0f} GB below --min-avail-gb')
    cap = int(args.memory_gb * (1024 ** 3))
    resource.setrlimit(resource.RLIMIT_AS, (cap, cap))
    for var in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        os.environ[var] = str(args.threads)
    return dict(mem_available_gb=round(avail, 1), rlimit_as_gb=args.memory_gb, threads=args.threads,
                loadavg=open('/proc/loadavg').read().split()[:3])

def mixture_tail_imhof(Q, lam):
    from scipy.integrate import quad
    lam = C.np.asarray(lam, dtype=float)
    lam = lam[lam > lam.max() * 1e-12] if lam.max() > 0 else lam

    def integrand(u):
        theta = 0.5 * C.np.sum(C.np.arctan(lam * u)) - 0.5 * Q * u
        rho = C.np.exp(0.25 * C.np.sum(C.np.log1p((lam * u) ** 2)))
        return math.sin(theta) / (u * rho)

    scale = 1.0 / max(lam.sum(), 1e-12)
    total, err = 0.0, 0.0
    for a, b in ((1e-10, scale), (scale, 20 * scale), (20 * scale, 2000 * scale)):
        v, e = quad(integrand, a, b, limit=400)
        total += v
        err += abs(e)
    p = 0.5 - total / math.pi
    if not (0.0 < p < 1.0) or err > 1e-4:
        raise ValueError('imhof integration did not converge')
    return p

def mixture_tail_liu(Q, lam):
    from scipy.stats import ncx2, chi2
    lam = C.np.asarray(lam, dtype=float)
    c = [float((lam ** k).sum()) for k in (1, 2, 3, 4)]
    if c[1] <= 0:
        return 1.0
    s1 = c[2] / c[1] ** 1.5
    s2 = c[3] / c[1] ** 2
    sigma_q = math.sqrt(2 * c[1])
    if s1 ** 2 > s2:
        a = 1 / (s1 - math.sqrt(s1 ** 2 - s2))
        delta = s1 * a ** 3 - a ** 2
        dof = a ** 2 - 2 * delta
        mu_x, sigma_x = dof + delta, math.sqrt(2) * a
        z = (Q - c[0]) / sigma_q * sigma_x + mu_x
        return float(ncx2.sf(max(z, 0.0), dof, delta))
    dof = 1 / s2
    mu_x, sigma_x = dof, math.sqrt(2 * dof)
    z = (Q - c[0]) / sigma_q * sigma_x + mu_x
    return float(chi2.sf(max(z, 0.0), dof))

def skat_p(Q, lam, mode, counters):
    if mode == 'imhof':
        try:
            return mixture_tail_imhof(Q, lam), 'imhof'
        except Exception:
            counters['imhof_fallback'] += 1
    return mixture_tail_liu(Q, lam), 'liu'

def acat(pvals):
    np = C.np
    p = np.asarray(pvals, dtype=np.longdouble)
    tiny = np.finfo(float).tiny
    p = np.clip(p, tiny, 1 - np.finfo(float).eps)
    T = np.mean(np.tan((np.longdouble(0.5) - p) * np.longdouble(np.pi)))
    out = float(np.longdouble(0.5) - np.arctan(T) / np.longdouble(np.pi))
    return min(max(out, tiny), 1.0), float(T)

def main():
    args = parser().parse_args()
    shuffle = None if args.resid_shuffle_seed == 'none' else int(args.resid_shuffle_seed)
    C.require(shuffle is not None or args.allow_real_residuals,
              'real residuals require --allow-real-residuals; the pre-registration is not sealed')
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    g = guard(args)
    C.load_libraries(args.threads)
    np = C.np
    from scipy.stats import chi2

    chroms = None if args.chrom == 'all' else [c.strip() for c in args.chrom.split(',') if c.strip()]
    traits = [t.strip() for t in args.traits.split(',') if t.strip()]
    betas = []
    for spec in args.beta_weights.split(';'):
        a, b = spec.split(',')
        betas.append((float(a), float(b)))
    C.require(bool(betas) and bool(traits), 'empty --beta-weights/--traits')
    cadd_dir = None if args.cadd_dir == 'none' else args.cadd_dir
    strata = None if args.perm_strata == 'none' else args.perm_strata

    config = dict(arguments=vars(args), code_sha256=C.sha_file(__file__),
                  common_sha256=C.sha_file(Path(__file__).with_name('l3_common.py')))
    C.atomic_json(out / f'{args.tag}.RUN.json', dict(**config, guard=g))

    cache = C.BandCache(args.root, cadd_dir, args.cadd_source, chroms, args.t1_partial_policy)
    ref_rows, M = C.ecdf_reference(cache, args.ecdf_reference)
    annots = ['__flat__'] + C.STAAR_ANNOTATIONS
    pi = {}
    ecdf_meta = {}
    for k in annots:
        if k == '__flat__':
            pi[k] = np.ones(len(cache.keys))
            ecdf_meta[k] = dict(definition='A = 1 (paper k=0 arm)')
        else:
            pi[k], ecdf_meta[k] = C.ecdf_avg_rank(cache.col(k), ref_rows, M)
    maf = cache.col('maf')
    C.require(np.isfinite(maf).all() and (maf > 0).all() and (maf <= .5).all(), 'invalid MAF for Beta weights')
    C.log(f'STAAR-O: rows={len(cache.keys)} genes={cache.ng} annotations={len(annots)} '
          f'beta_weights={betas} ECDF_M={M}')

    counters = dict(imhof_fallback=0, skipped_large=0, skipped_singular=0)
    rows, resid, perms, resid_samples = [], None, None, None
    ds_fp = {}
    for chrom in cache.use_chr:
        ds = C.Dosage(args.root, chrom, samples=resid_samples)
        ds_fp.update(ds.fingerprint)
        if resid is None:
            resid_samples = ds.samples
            resid = C.Residuals(args.resid_dir, traits, resid_samples, args.min_n, shuffle, strata)
            if args.staar_null == 'permutation':
                perms = resid.permutations(args.perm_b, args.perm_seed)
            C.log(f'residuals shuffled={resid.shuffled} blocks={resid.meta["n_blocks"]}')
        gsel = [gi for gi in range(cache.ng) if cache.gchr[cache.gl[gi]] == chrom]
        C.log(f'chr{chrom}: genes={len(gsel)}')
        for gi in gsel:
            pp = cache.pairs_of[gi]
            keys = cache.keys[pp]
            uk, ui = np.unique(keys, return_index=True)
            pp = pp[np.sort(ui)]
            if len(pp) > args.max_variants_per_gene:
                counters['skipped_large'] += 1
                rows.append([chrom, cache.gl[gi], 'ALL', len(pp), '', '', '', '', 'skipped_large'])
                continue
            vloc = ds.locate(cache.keys[pp])
            Gs = ds.M[vloc]
            m = Gs.shape[0]
            colsum_all = np.asarray(Gs.sum(axis=1)).ravel()
            for it, trait in enumerate(traits):
                obs = np.flatnonzero(resid.mask[it])
                Go = Gs[:, obs]
                n = len(obs)
                r = resid.full[it, obs].astype(np.float64)
                r = r - r.mean()
                sigma2 = float((r * r).sum() / n)
                cs = np.asarray(Go.sum(axis=1)).ravel()
                GG = np.asarray((Go @ Go.T).todense(), dtype=np.float64)
                Sig = sigma2 * (GG - np.outer(cs, cs) / n)
                U = np.asarray(Go @ r, dtype=np.float64)
                pset, plabels = [], []
                for (ba, bb) in betas:
                    wb = np.exp((ba - 1) * np.log(maf[pp]) + (bb - 1) * np.log1p(-maf[pp])
                                - math.lgamma(ba) - math.lgamma(bb) + math.lgamma(ba + bb))
                    for k in annots:
                        wt = pi[k][pp] * wb
                        vb = float(wt @ Sig @ wt)
                        if vb <= 0:
                            counters['skipped_singular'] += 1
                            continue
                        qb = float(wt @ U) ** 2 / vb
                        pb = float(chi2.sf(qb, 1))
                        A = Sig * wt[:, None] * wt[None, :]
                        qs = float(np.sum((wt * U) ** 2))
                        lam = np.linalg.eigvalsh((A + A.T) / 2)
                        lam = lam[lam > 0]
                        ps = 1.0 if len(lam) == 0 else skat_p(qs, lam, args.skat_tail, counters)[0]
                        pset.append(pb); plabels.append(f'burden|{k}|B{ba:g},{bb:g}')
                        pset.append(ps); plabels.append(f'skat|{k}|B{ba:g},{bb:g}')
                if not pset:
                    rows.append([chrom, cache.gl[gi], trait, m, '', '', '', '', 'no_testable_weight'])
                    continue
                keep = [i for i, lab in enumerate(plabels)
                        if (args.acat_set == 'annot-x-weight-x-test')
                        or (args.acat_set == 'annot-x-weight' and lab.startswith('skat|'))
                        or (args.acat_set == 'annot-only-burden' and lab.startswith('burden|')
                            and lab.endswith('B1,1'))]
                C.require(bool(keep), 'the requested --acat-set selected no p-values')
                p_o, tstat = acat([pset[i] for i in keep])
                rows.append([chrom, cache.gl[gi], trait, m, len(keep),
                             float(min(pset[i] for i in keep)), float(p_o), float(tstat), 'ok'])
                if len(rows) % 200 == 0:
                    C.log(f'  {len(rows)} gene/trait records')
        del ds
        C.log(f'chr{chrom} done; records={len(rows)}')

    C.require(bool(rows), 'no STAAR-O records produced')
    header = ['chr', 'gene', 'trait', 'n_variants', 'n_pvalues_combined', 'p_min_component',
              'p_STAAR_O', 'acat_T', 'status']
    C.atomic_tsv(out / f'{args.tag}.staaro.tsv', header, rows)
    ok = [r for r in rows if r[8] == 'ok']
    summary = dict(config=config, guard=g, annotations=annots, ecdf=ecdf_meta, ecdf_M=int(M),
                   beta_weights=[list(b) for b in betas], acat_set=args.acat_set,
                   staar_null=args.staar_null, skat_tail=args.skat_tail, counters=counters,
                   residuals=resid.meta, n_records=len(rows), n_ok=len(ok),
                   p_summary=(dict(min=float(min(r[6] for r in ok)),
                                   median=float(C.np.median([r[6] for r in ok])),
                                   frac_below_0p05=float(C.np.mean([r[6] < .05 for r in ok])))
                              if ok else {}),
                   input_fingerprints=dict(**cache.fingerprints, **ds_fp),
                   residual_mode='SHUFFLED (mechanical smoke)' if resid.shuffled
                                 else 'REAL residuals — sealed pre-registration required',
                   fidelity_notes=[
                       'pi_hat and MAF weights multiply the score vector; the Cauchy/ACAT transform '
                       'is applied to p-values only.',
                       'Null covariance uses sigma^2 * Gc\' Gc on covariate-adjusted residuals; the '
                       'covariate projection P is not recoverable from the residual files.',
                       'ACAT-V is not part of this arm — the task specifies Burden and SKAT only.'])
    C.atomic_json(out / f'{args.tag}.SUMMARY.json', summary)
    C.atomic_json(out / f'{args.tag}.done',
                  dict(status='complete', tag=args.tag, records=len(rows),
                       table_sha256=C.sha_file(out / f'{args.tag}.staaro.tsv'),
                       residual_shuffled=bool(resid.shuffled)))
    C.log('L3_STAARO_DONE')

if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        msg = str(error) if isinstance(error, RuntimeError) and str(error).startswith('GATE FAIL:') \
            else type(error).__name__
        print(f'[l3_staaro] {msg}', flush=True)
        raise SystemExit(1) from None
