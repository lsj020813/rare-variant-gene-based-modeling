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
import csv
import json
import math
from pathlib import Path

import l3_common as C

def parser():
    p = argparse.ArgumentParser(description='layer-3 sealed metric table')
    p.add_argument('--gene-trait', required=True, help='l3_burden <tag>.gene_trait.tsv')
    p.add_argument('--burden-summary', required=True, help='l3_burden <tag>.SUMMARY.json')
    p.add_argument('--staaro', required=True, help="l3_staaro <tag>.staaro.tsv, or 'none'")
    p.add_argument('--truth-dir', required=True, help="SAIGE-GENE result directory, or 'none'")
    p.add_argument('--truth-group', required=True, help="Group value to keep, e.g. 'all'")
    p.add_argument('--truth-maxmaf', required=True, help="max_MAF value to keep, e.g. '0.01'")
    p.add_argument('--truth-alpha', required=True, type=float, help='truth significance threshold')
    p.add_argument('--reference-arm', required=True, help="denominator arm for (na), e.g. 'flat'")
    p.add_argument('--alpha-grid', required=True, help='comma list of thresholds for (ga) and (ra)')
    p.add_argument('--topk-grid', required=True, help='comma list of k for (da)')
    p.add_argument('--p-column', required=True, choices=['p', 'p_adj'], help='which p drives (ga)/(da)/(ra)')
    p.add_argument('--out', required=True)
    p.add_argument('--tag', required=True)
    return p

def read_tsv(path):
    with open(path, newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))

def acat(pvals):
    np = C.np
    p = np.asarray(pvals, dtype=np.longdouble)
    tiny = np.finfo(float).tiny
    p = np.clip(p, tiny, 1 - np.finfo(float).eps)
    T = np.mean(np.tan((np.longdouble(0.5) - p) * np.longdouble(np.pi)))
    return float(min(max(float(np.longdouble(0.5) - np.arctan(T) / np.longdouble(np.pi)), tiny), 1.0))

def load_truth(truth_dir, traits, group, maxmaf, alpha, chroms):
    truth, meta = {}, {}
    for trait in traits:
        sig, tested, missing = set(), set(), []
        for ch in chroms:
            f = Path(truth_dir) / f'{trait}.chr{ch}'
            if not (f.is_file() and Path(str(f) + '.done').is_file()):
                missing.append(ch)
                continue
            for row in read_tsv(f):
                if row.get('Group') != group or row.get('max_MAF') != maxmaf:
                    continue
                gene = row['Region'].split('.')[0]
                tested.add(gene)
                try:
                    if float(row['Pvalue']) < alpha:
                        sig.add(gene)
                except (TypeError, ValueError):
                    continue
        truth[trait] = dict(sig=sig, tested=tested)
        meta[trait] = dict(n_truth_tested=len(tested), n_truth_significant=len(sig),
                           chromosomes_missing=missing)
    return truth, meta

def main():
    args = parser().parse_args()
    C.load_libraries(1)
    np = C.np
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    table = read_tsv(args.gene_trait)
    C.require(bool(table), 'empty gene_trait table')
    with open(args.burden_summary) as fh:
        bsum = json.load(fh)
    alphas = [float(x) for x in args.alpha_grid.split(',') if x.strip()]
    topks = [int(x) for x in args.topk_grid.split(',') if x.strip()]
    arms = sorted({r['arm'] for r in table})
    traits = sorted({r['trait'] for r in table})
    chroms = sorted({r['chr'] for r in table}, key=int)
    C.require(args.reference_arm in arms, 'reference arm absent from the table')

    cell = {}
    for r in table:
        cell.setdefault((r['arm'], r['trait']), {})[r['gene']] = (
            float(r['T']), float(r[args.p_column]), int(r['n_unique_variants']),
            int(r['n_phi_feature_complete']))
    for arm in arms:
        genes = sorted(set.intersection(*[set(cell[(arm, t)]) for t in traits])) if traits else []
        comb = {}
        for gene in genes:
            ps = [cell[(arm, t)][gene][1] for t in traits]
            ts = [abs(cell[(arm, t)][gene][0]) for t in traits]
            comb[gene] = (float(np.mean(ts)), acat(ps),
                          cell[(arm, traits[0])][gene][2], cell[(arm, traits[0])][gene][3])
        cell[(arm, 'COMBINED')] = comb
    report_traits = traits + ['COMBINED']

    truth, tmeta = ({}, {})
    if args.truth_dir != 'none':
        truth, tmeta = load_truth(args.truth_dir, traits, args.truth_group, args.truth_maxmaf,
                                  args.truth_alpha, chroms)

    staaro = {}
    if args.staaro != 'none':
        for r in read_tsv(args.staaro):
            if r['status'] == 'ok':
                staaro.setdefault(r['trait'], {})[r['gene']] = float(r['p_STAAR_O'])

    rows = []
    for arm in arms:
        for trait in report_traits:
            d = cell[(arm, trait)]
            genes = sorted(d)
            if not genes:
                continue
            T = np.array([d[g][0] for g in genes])
            p = np.array([d[g][1] for g in genes])
            ref = cell[(args.reference_arm, trait)]
            Tref = np.array([ref[g][0] for g in genes if g in ref])
            pref = np.array([ref[g][1] for g in genes if g in ref])

            mean_abs = float(np.abs(T).mean())
            mean_abs_ref = float(np.abs(Tref).mean()) if len(Tref) else float('nan')
            ratio = mean_abs / mean_abs_ref if mean_abs_ref > 0 else None
            pz = (bsum.get('primary_metric', {}) or {}).get(f'{arm}|{trait}', {})
            perm_z = pz.get('perm_z')

            counts = {}
            for a in alphas:
                n_arm = int((p <= a).sum())
                n_ref = int((pref <= a).sum()) if len(pref) else 0
                counts[f'{a:g}'] = dict(arm=n_arm, reference=n_ref,
                                        verdict=('unmeasurable' if (n_arm == 0 and n_ref == 0)
                                                 else ('gain' if n_arm > n_ref
                                                       else ('loss' if n_arm < n_ref else 'tie'))))

            order = np.argsort(p, kind='stable')
            recall = {}
            if trait != 'COMBINED' and trait in truth:
                tset = truth[trait]['sig'] & set(genes)
                for k in topks:
                    top = {genes[i] for i in order[:k]}
                    recall[str(k)] = dict(
                        k=k, truth_significant_in_scope=len(tset),
                        hit=len(top & tset),
                        recall=(len(top & tset) / len(tset)) if tset else None,
                        note=None if tset else 'no truth-significant gene in scope')
            rows.append(dict(
                arm=arm, trait=trait, n_genes=len(genes),
                metric_na_mean_abs_T=mean_abs,
                metric_na_mean_abs_T_reference=mean_abs_ref,
                metric_na_ratio=ratio, metric_na_perm_z=perm_z,
                metric_na_mde_reference=MDE_REFERENCE,
                metric_ga_counts=counts,
                metric_da_topk_recall=recall,
                p_min=float(p.min()), T_abs_max=float(np.abs(T).max()),
                variants_per_gene_median=float(np.median([d[g][2] for g in genes])),
                phi_feature_complete_median=float(np.median([d[g][3] for g in genes])),
                staaro_p_min=(float(min(staaro[trait].values()))
                              if trait in staaro and staaro[trait] else None)))

    curve = []
    for r in rows:
        for a, v in r['metric_ga_counts'].items():
            curve.append([r['arm'], r['trait'], a, v['arm'], v['reference'], v['verdict']])
    C.atomic_tsv(out / f'{args.tag}.metric_ra_threshold_curve.tsv',
                 ['arm', 'trait', 'threshold', 'n_significant_arm', 'n_significant_reference', 'verdict'],
                 curve)

    flat_rows = []
    for r in rows:
        flat_rows.append([r['arm'], r['trait'], r['n_genes'],
                          r['metric_na_ratio'], r['metric_na_perm_z'], r['metric_na_mean_abs_T'],
                          r['metric_na_mean_abs_T_reference'],
                          json.dumps({k: v['arm'] for k, v in r['metric_ga_counts'].items()}),
                          json.dumps({k: v['reference'] for k, v in r['metric_ga_counts'].items()}),
                          json.dumps({k: v['verdict'] for k, v in r['metric_ga_counts'].items()}),
                          json.dumps({k: v.get('recall') for k, v in r['metric_da_topk_recall'].items()}),
                          json.dumps({k: v.get('hit') for k, v in r['metric_da_topk_recall'].items()}),
                          r['p_min'], r['staaro_p_min'],
                          r['variants_per_gene_median'], r['phi_feature_complete_median']])
    C.atomic_tsv(out / f'{args.tag}.metrics_four.tsv',
                 ['arm', 'trait', 'n_genes', 'na_ratio_vs_reference', 'na_perm_z',
                  'na_mean_abs_T', 'na_mean_abs_T_reference', 'ga_n_sig_arm', 'ga_n_sig_reference',
                  'ga_verdict', 'da_topk_recall', 'da_topk_hits', 'p_min', 'staaro_p_min',
                  'median_variants_per_gene', 'median_phi_feature_complete'],
                 flat_rows)
    C.atomic_json(out / f'{args.tag}.metrics.json',
                  dict(config=vars(args), reference_arm=args.reference_arm,
                       residual_mode=bsum.get('residual_mode'),
                       permutation=bsum.get('permutation_meta'), truth=tmeta, rows=rows,
                       mde_reference=MDE_REFERENCE,
                       reading_rules=[
                           '(ga) 0 vs 0 is reported as unmeasurable, never as no gain.',
                           '(na) is the sealed primary metric; (ra) is supporting evidence.',
                           'COMBINED trait = ACAT of the per-trait p-values within a gene.']))
    C.atomic_json(out / f'{args.tag}.done', dict(status='complete', rows=len(flat_rows)))
    C.log('L3_METRICS_DONE')

if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        msg = str(error) if isinstance(error, RuntimeError) and str(error).startswith('GATE FAIL:') \
            else type(error).__name__
        print(f'[l3_metrics] {msg}', flush=True)
        raise SystemExit(1) from None
