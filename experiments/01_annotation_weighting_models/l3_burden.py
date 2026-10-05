#!/usr/bin/env python3
import os as _cfg_os
import math as _cfg_math

def _config_number(name, cast=int, positive=False):
    value = cast(_cfg_os.environ[name])
    if not _cfg_math.isfinite(value) or (positive and value <= 0):
        raise ValueError(name + " must be finite" + (" and positive" if positive else ""))
    return value

import os as _os
MDE_REFERENCE = _config_number("MDE_REFERENCE", float, True)
import argparse
import json
import os
from pathlib import Path
import resource
import sys
import time

import l3_common as C

def parser():
    p = argparse.ArgumentParser(description='layer-3 gene burden test (C2)')
    p.add_argument('--root', required=True, help='band reference root (holds annot/cache, annot/ds, groupfiles_bwg)')
    p.add_argument('--cadd-dir', required=True, help="directory holding the cadd source, or 'none'")
    p.add_argument('--cadd-source', required=True, choices=['features', 'extract'],
                   help='features = chrN.features.tsv.gz:cadd_phred (22/22); extract = chrN.cadd.tsv (12/22)')
    p.add_argument('--chrom', required=True, help="comma list of chromosomes, or 'all'")
    p.add_argument('--traits', required=True, help='comma list of trait codes')
    p.add_argument('--resid-dir', required=True, help='build_resid_v8.py output directory')
    p.add_argument('--arms', required=True, help='comma list: learned,flat,cadd_rank,pibar')
    p.add_argument('--phi-json', required=True, help="frozen phi export, or 'none'")
    p.add_argument('--v10-script', required=True, help="l1_train_v10.py path, or 'none'")
    p.add_argument('--pair-weight', required=True, choices=['none', 'inverse-n-genes'],
                   help='pw_v; v8 used inverse-n-genes')
    p.add_argument('--ecdf-reference', required=True, choices=['band-all', 'band-chrom', 'band-cadd-available'])
    p.add_argument('--t1-partial-policy', required=True, choices=['blank', 'fail'])
    p.add_argument('--perm-b', required=True, type=int, help='number of permutations B')
    p.add_argument('--perm-seed', required=True, type=int)
    p.add_argument('--perm-block', required=True, type=int, help='permutations per GEMM chunk')
    p.add_argument('--pmode', required=True, choices=['empirical', 'moment'])
    p.add_argument('--alpha', required=True, type=float)
    p.add_argument('--multiple', required=True,
                   choices=['none', 'bonferroni-within-trait', 'bonferroni-all', 'bh-within-trait', 'bh-all'])
    p.add_argument('--min-n', required=True, type=int)
    p.add_argument('--gene-batch', required=True, type=int)
    p.add_argument('--threads', required=True, type=int)
    p.add_argument('--memory-gb', required=True, type=float, help='RLIMIT_AS ceiling for this process')
    p.add_argument('--min-avail-gb', required=True, type=float, help='refuse to start below this MemAvailable')
    p.add_argument('--resid-shuffle-seed', required=True, help="int seed, or 'none' for real residuals")
    p.add_argument('--perm-strata', required=True, help="IID/stratum TSV, or 'none'")
    p.add_argument('--allow-real-residuals', action='store_true',
                   help='required to run with --resid-shuffle-seed none (sealed pre-registration only)')
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
    ncore = sum(1 for l in open('/proc/cpuinfo') if l.startswith('processor'))
    return dict(mem_available_gb=round(avail, 1), rlimit_as_gb=args.memory_gb,
                threads=args.threads, machine_cores=ncore, loadavg=open('/proc/loadavg').read().split()[:3])

def bh(pvals):
    order = C.np.argsort(pvals)
    m = len(pvals)
    ranked = pvals[order] * m / C.np.arange(1, m + 1)
    adj = C.np.minimum.accumulate(ranked[::-1])[::-1]
    out = C.np.empty(m)
    out[order] = C.np.minimum(adj, 1.0)
    return out

def adjust(mode, pvals, trait_index, alpha):
    p = C.np.asarray(pvals, dtype=float)
    out = p.copy()
    if mode == 'none':
        return out
    if mode == 'bonferroni-all':
        return C.np.minimum(p * len(p), 1.0)
    if mode == 'bh-all':
        return bh(p)
    for t in C.np.unique(trait_index):
        sel = trait_index == t
        out[sel] = C.np.minimum(p[sel] * int(sel.sum()), 1.0) if mode == 'bonferroni-within-trait' else bh(p[sel])
    return out

def main():
    args = parser().parse_args()
    C.require(args.perm_b >= 1 and args.perm_block >= 1 and args.gene_batch >= 1, 'perm/batch sizes must be positive')
    shuffle = None if args.resid_shuffle_seed == 'none' else int(args.resid_shuffle_seed)
    C.require(shuffle is not None or args.allow_real_residuals,
              'real residuals require --allow-real-residuals; the pre-registration is not sealed')
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    g = guard(args)
    C.load_libraries(args.threads)
    np = C.np
    from scipy.stats import norm

    chroms = None if args.chrom == 'all' else [c.strip() for c in args.chrom.split(',') if c.strip()]
    traits = [t.strip() for t in args.traits.split(',') if t.strip()]
    arms = [a.strip() for a in args.arms.split(',') if a.strip()]
    C.require(bool(traits) and bool(arms), 'empty --traits/--arms')
    cadd_dir = None if args.cadd_dir == 'none' else args.cadd_dir
    v10 = None if args.v10_script == 'none' else C.import_v10(args.v10_script)
    phi_json = None if args.phi_json == 'none' else args.phi_json
    strata = None if args.perm_strata == 'none' else args.perm_strata

    config = dict(arguments=vars(args), code_sha256=C.sha_file(__file__),
                  common_sha256=C.sha_file(Path(__file__).with_name('l3_common.py')))
    C.atomic_json(out / f'{args.tag}.RUN.json', dict(**config, guard=g))

    cache = C.BandCache(args.root, cadd_dir, args.cadd_source, chroms, args.t1_partial_policy)
    C.log(f'cache rows={len(cache.keys)} genes={cache.ng} chroms={",".join(cache.use_chr)} '
          f'cadd_chroms={len(cache.cadd_available)} t1_partial_blanked={cache.t1_partial_blanked}')

    pw = cache.pair_weight(args.pair_weight).astype(np.float64)
    phis, phi_meta = {}, {}
    for arm in arms:
        phis[arm], phi_meta[arm] = C.phi_arm(cache, arm, v10, phi_json, args.ecdf_reference)
        v = phis[arm]
        phi_meta[arm]['summary'] = dict(n=int(len(v)), min=float(v.min()), max=float(v.max()),
                                        mean=float(v.mean()), sd=float(v.std()))
        C.log(f'arm {arm}: phi mean={v.mean():.6g} sd={v.std():.6g} range=[{v.min():.6g},{v.max():.6g}]')

    t1_ok = cache.col('t1_na') == 0
    cadd_ok = np.isfinite(cache.col('cadd'))
    coverage = dict(band_pairs=int(len(cache.keys)),
                    band_unique_variants=int(len(cache.unique_rows)),
                    pairs_t1_observed=int(t1_ok.sum()), pairs_cadd_observed=int(cadd_ok.sum()),
                    pairs_phi_feature_block_complete=int((t1_ok & cadd_ok).sum()),
                    unique_t1_observed=int(t1_ok[cache.unique_rows].sum()),
                    unique_cadd_observed=int(cadd_ok[cache.unique_rows].sum()))

    rows, gene_meta = [], []
    resid = None
    perms = None
    ds_fp = {}
    abs_obs, abs_perm, n_gene_counted = {}, {}, 0
    for chrom in cache.use_chr:
        ds = C.Dosage(args.root, chrom, samples=None if resid is None else resid_samples)
        ds_fp.update(ds.fingerprint)
        if resid is None:
            resid_samples = ds.samples
            resid = C.Residuals(args.resid_dir, traits, resid_samples, args.min_n, shuffle, strata)
            R = resid.R.astype(np.float32)
            nt = len(traits)
            perms = resid.permutations(args.perm_b, args.perm_seed)
            C.log(f'residuals loaded shuffled={resid.shuffled} blocks={resid.meta["n_blocks"]} '
                  f'movable={resid.meta["movable_observed_samples"]} B={args.perm_b}')
        else:
            C.require(ds.samples == resid_samples, 'chromosome sample order mismatch')
        gsel = [gi for gi in range(cache.ng) if cache.gchr[cache.gl[gi]] == chrom]
        C.log(f'chr{chrom}: dosage variants={ds.n_variants} nnz={ds.nnz} genes={len(gsel)}')
        for lo in range(0, len(gsel), args.gene_batch):
            block = gsel[lo:lo + args.gene_batch]
            vlocs = [ds.locate(cache.keys[cache.pairs_of[gi]]) for gi in block]
            for arm in arms:
                w = phis[arm] * pw
                S = np.empty((len(block), ds.ns), dtype=np.float64)
                for bi, gi in enumerate(block):
                    S[bi] = C.burden(ds, vlocs[bi], w[cache.pairs_of[gi]].astype(np.float32))
                Z, sd = C.zrow(S)
                C.require(np.isfinite(Z).all(), 'nonfinite Z')
                T = Z @ R.T.astype(np.float64)
                ge = np.zeros_like(T, dtype=np.int64)
                psum = np.zeros_like(T)
                psq = np.zeros_like(T)
                abs_obs.setdefault(arm, np.zeros(nt))
                abs_perm.setdefault(arm, np.zeros((args.perm_b, nt)))
                abs_obs[arm] += np.abs(T).sum(axis=0)
                for pl in range(0, args.perm_b, args.perm_block):
                    chunk = perms[pl:pl + args.perm_block]
                    Rp = np.concatenate([R[:, pm] for pm in chunk], axis=0).astype(np.float64)
                    Tp = Z @ Rp.T
                    Tp = Tp.reshape(len(block), len(chunk), nt)
                    ge += (np.abs(Tp) >= np.abs(T)[:, None, :]).sum(axis=1)
                    psum += Tp.sum(axis=1)
                    psq += (Tp * Tp).sum(axis=1)
                    abs_perm[arm][pl:pl + len(chunk)] += np.abs(Tp).sum(axis=0)
                    del Rp, Tp
                pmean = psum / args.perm_b
                pvar = np.maximum(psq / args.perm_b - pmean ** 2, 0.0)
                psd = np.sqrt(pvar)
                if args.pmode == 'empirical':
                    pv = (1.0 + ge) / (args.perm_b + 1.0)
                else:
                    C.require((psd > 0).all(), 'zero permutation sd under --pmode moment')
                    pv = 2.0 * norm.sf(np.abs((T - pmean) / psd))
                for bi, gi in enumerate(block):
                    pp = cache.pairs_of[gi]
                    nv, nv_uniq = len(pp), len(set(cache.keys[pp]))
                    n_feat = int((t1_ok[pp] & cadd_ok[pp]).sum())
                    for it, trait in enumerate(traits):
                        rows.append([chrom, cache.gl[gi], arm, trait, nv, nv_uniq, n_feat,
                                     float(sd[bi, 0]), float(T[bi, it]), float(pmean[bi, it]),
                                     float(psd[bi, it]), int(ge[bi, it]), float(pv[bi, it])])
                    if arm == arms[0]:
                        gene_meta.append(dict(chr=chrom, gene=cache.gl[gi], n_pairs=nv,
                                              n_unique_variants=nv_uniq, n_phi_feature_complete=n_feat))
                        n_gene_counted += 1
                del S, Z, T, ge, psum, psq
        del ds
        C.log(f'chr{chrom} done; rows so far={len(rows)}')

    C.require(bool(rows), 'no gene rows produced')
    pv = np.array([r[12] for r in rows])
    tix = np.array([traits.index(r[3]) * 100 + arms.index(r[2]) for r in rows])
    padj = adjust(args.multiple, pv, tix, args.alpha)
    header = ['chr', 'gene', 'arm', 'trait', 'n_pairs', 'n_unique_variants', 'n_phi_feature_complete',
              'S_sd', 'T', 'T_perm_mean', 'T_perm_sd', 'n_perm_ge', 'p', 'p_adj', 'significant']
    table = [r + [float(padj[i]), int(padj[i] <= args.alpha)] for i, r in enumerate(rows)]
    C.atomic_tsv(out / f'{args.tag}.gene_trait.tsv', header, table)

    primary = {}
    if 'flat' in abs_obs:
        ng = max(n_gene_counted, 1)
        for arm in arms:
            for it, trait in enumerate(traits):
                num_o, den_o = abs_obs[arm][it], abs_obs['flat'][it]
                ratio = float(num_o / den_o) if den_o > 0 else None
                nullr = abs_perm[arm][:, it] / np.maximum(abs_perm['flat'][:, it], np.finfo(float).tiny)
                mu, sd = float(nullr.mean()), float(nullr.std())
                primary[f'{arm}|{trait}'] = dict(
                    metric='mean_abs_T_ratio_vs_flat', genes=int(ng),
                    mean_abs_T=float(num_o / ng), mean_abs_T_flat=float(den_o / ng),
                    ratio=ratio, perm_ratio_mean=mu, perm_ratio_sd=sd,
                    perm_z=(float((ratio - mu) / sd) if (ratio is not None and sd > 0) else None),
                    B=int(args.perm_b),
                    note='L3-1 sealed primary metric; MDE reference MDE_REFERENCE')
    else:
        primary = dict(unavailable='the flat arm must be included to form the L3-1 ratio')

    summary = dict(config=config, guard=g, coverage=coverage, primary_metric=primary,
                   residuals=resid.meta, permutation_meta=dict(B=args.perm_b, seed=args.perm_seed,
                                                               pmode=args.pmode, block=args.perm_block),
                   arms=phi_meta, chromosomes=cache.use_chr, n_genes=cache.ng,
                   gene_variant_distribution=gene_variant_stats(gene_meta),
                   per_arm_trait=per_arm_trait_stats(rows, arms, traits, args.alpha, padj),
                   input_fingerprints=dict(**cache.fingerprints, **ds_fp),
                   residual_mode='SHUFFLED (mechanical smoke)' if resid.shuffled
                                 else 'REAL residuals — sealed pre-registration required')
    C.atomic_json(out / f'{args.tag}.SUMMARY.json', summary)
    C.atomic_json(out / f'{args.tag}.done',
                  dict(status='complete', tag=args.tag, rows=len(table),
                       table_sha256=C.sha_file(out / f'{args.tag}.gene_trait.tsv'),
                       summary_sha256=C.sha_file(out / f'{args.tag}.SUMMARY.json'),
                       residual_shuffled=bool(resid.shuffled)))
    C.log('L3_BURDEN_DONE')

def gene_variant_stats(gene_meta):
    np = C.np
    if not gene_meta:
        return {}
    n = np.array([g['n_pairs'] for g in gene_meta])
    u = np.array([g['n_unique_variants'] for g in gene_meta])
    f = np.array([g['n_phi_feature_complete'] for g in gene_meta])
    q = [0, 5, 25, 50, 75, 95, 100]
    return dict(genes=len(gene_meta),
                pairs_per_gene={f'p{x}': float(np.percentile(n, x)) for x in q},
                unique_variants_per_gene={f'p{x}': float(np.percentile(u, x)) for x in q},
                phi_feature_complete_per_gene={f'p{x}': float(np.percentile(f, x)) for x in q},
                genes_with_zero_feature_complete=int((f == 0).sum()),
                total_pairs=int(n.sum()), total_unique=int(u.sum()))

def per_arm_trait_stats(rows, arms, traits, alpha, padj):
    np = C.np
    T = np.array([r[8] for r in rows])
    p = np.array([r[12] for r in rows])
    out = {}
    for arm in arms:
        for trait in traits:
            sel = np.array([(r[2] == arm and r[3] == trait) for r in rows])
            if not sel.any():
                continue
            t, pp, pa = T[sel], p[sel], padj[sel]
            out[f'{arm}|{trait}'] = dict(
                n_genes=int(sel.sum()), T_mean=float(t.mean()), T_sd=float(t.std()),
                T_median=float(np.median(t)), T_abs_max=float(np.abs(t).max()),
                p_min=float(pp.min()), frac_p_below_0p05=float((pp < .05).mean()),
                frac_p_below_0p5=float((pp < .5).mean()), n_significant=int((pa <= alpha).sum()))
    return out

if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        msg = str(error) if isinstance(error, RuntimeError) and str(error).startswith('GATE FAIL:') \
            else type(error).__name__
        print(f'[l3_burden] {msg}', flush=True)
        raise SystemExit(1) from None
