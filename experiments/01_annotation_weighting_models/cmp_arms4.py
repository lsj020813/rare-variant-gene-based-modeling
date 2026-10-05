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


import argparse, csv, json, math
from pathlib import Path
import sys
sys.path.insert(0, _config_path('${PROJECT_ROOT}/work/run_l3'))
import l3_common as C

ARMS = ['none', 'flat', 'phi', 'phibeta']

def load(path):
    d = {}
    with open(path, newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            try:
                d[r['Region']] = dict(
                    P=float(r['Pvalue']),
                    PB=float(r['Pvalue_Burden']) if r.get('Pvalue_Burden') not in (None, '', 'NA') else None,
                    PS=float(r['Pvalue_SKAT']) if r.get('Pvalue_SKAT') not in (None, '', 'NA') else None,
                    MAC=float(r['MAC']) if r.get('MAC') not in (None, '', 'NA') else None,
                    nrare=r.get('Number_rare'))
            except (TypeError, ValueError):
                continue
    return d

def pair_stats(np, spearmanr, kendalltau, A, B, genes, tiny):
    a = np.array([A[g]['P'] for g in genes])
    b = np.array([B[g]['P'] for g in genes])
    d = -np.log10(np.maximum(a, tiny)) + np.log10(np.maximum(b, tiny))
    sp, kt = spearmanr(a, b), kendalltau(a, b)
    return dict(spearman_rho=float(sp.statistic), spearman_p=float(sp.pvalue),
                kendall_tau=float(kt.statistic),
                dlog10P=dict(mean=float(d.mean()), sd=float(d.std()), median=float(np.median(d)),
                             q={f'p{q}': float(np.percentile(d, q)) for q in (1, 5, 25, 50, 75, 95, 99)},
                             frac_first_more_significant=float((d > 0).mean()),
                             frac_abs_gt_0p5=float((np.abs(d) > 0.5).mean()),
                             frac_abs_gt_1=float((np.abs(d) > 1).mean()),
                             max_gain=float(d.max()), max_loss=float(d.min())))

def main():
    ap = argparse.ArgumentParser()
    for a in ARMS:
        ap.add_argument(f'--{a}', required=True)
    ap.add_argument('--truth', required=True)
    ap.add_argument('--truth-alpha', required=True, type=float)
    ap.add_argument('--thresholds', required=True)
    ap.add_argument('--just-below', required=True, help='lo,hi -- the (d) band, e.g. 2.5e-6,1e-5')
    ap.add_argument('--par-record', required=True, help='JSON of arm -> xargs parallelism actually used')
    ap.add_argument('--out', required=True)
    ap.add_argument('--tag', required=True)
    args = ap.parse_args()
    C.load_libraries(1)
    np = C.np
    from scipy.stats import spearmanr, kendalltau
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    tiny = np.finfo(float).tiny

    arms = {a: load(getattr(args, a)) for a in ARMS}
    truth = load(args.truth)
    genes = sorted(set.intersection(*[set(v) for v in arms.values()]) & set(truth))
    C.require(bool(genes), 'no genes shared by the four arms and the truth file')
    thr = [float(x) for x in args.thresholds.split(',')]
    lo, hi = (float(x) for x in args.just_below.split(','))

    gate = dict(genes_compared=len(genes),
                none_equals_truth=all(arms['none'][g]['P'] == truth[g]['P'] for g in genes),
                worst_abs_diff=max(abs(arms['none'][g]['P'] - truth[g]['P']) for g in genes))

    counts = {a: {f'{t:g}': int(sum(arms[a][g]['P'] <= t for g in genes)) for t in thr} for a in ARMS}
    counts['truth'] = {f'{t:g}': int(sum(truth[g]['P'] <= t for g in genes)) for t in thr}

    questions = dict(
        mechanism_phi_vs_flat=pair_stats(np, spearmanr, kendalltau, arms['phi'], arms['flat'], genes, tiny),
        paper_phibeta_vs_none=pair_stats(np, spearmanr, kendalltau, arms['phibeta'], arms['none'], genes, tiny),
        aux_phibeta_vs_phi=pair_stats(np, spearmanr, kendalltau, arms['phibeta'], arms['phi'], genes, tiny),
        aux_none_vs_flat=pair_stats(np, spearmanr, kendalltau, arms['none'], arms['flat'], genes, tiny))

    tsig = sorted((g for g in genes if truth[g]['P'] < args.truth_alpha), key=lambda g: truth[g]['P'])
    trows = []
    for g in tsig:
        r = dict(gene=g, truth_P=truth[g]['P'])
        for a in ARMS:
            r[f'{a}_P'] = arms[a][g]['P']
        r['dlog10P_phibeta_vs_none'] = (-math.log10(max(r['phibeta_P'], tiny))
                                        + math.log10(max(r['none_P'], tiny)))
        r['dlog10P_phi_vs_flat'] = (-math.log10(max(r['phi_P'], tiny))
                                    + math.log10(max(r['flat_P'], tiny)))
        r['phibeta_lower_than_none'] = bool(r['phibeta_P'] < r['none_P'])
        r['phi_lower_than_flat'] = bool(r['phi_P'] < r['flat_P'])
        trows.append(r)

    band = dict(lo=lo, hi=hi,
                definition=f'genes with {lo:g} < P <= {hi:g} in the reference arm')
    for a in ARMS:
        band[f'{a}_in_band'] = sorted((g for g in genes if lo < arms[a][g]['P'] <= hi),
                                      key=lambda g: arms[a][g]['P'])
    for tag, test, ref in (('phibeta_vs_none', 'phibeta', 'none'), ('phi_vs_flat', 'phi', 'flat')):
        crossed = sorted(g for g in genes if arms[ref][g]['P'] > lo >= arms[test][g]['P'])
        lost = sorted(g for g in genes if arms[test][g]['P'] > lo >= arms[ref][g]['P'])
        band[tag] = dict(crossed_into_significance=crossed, n_crossed=len(crossed),
                         lost_significance=lost, n_lost=len(lost),
                         band_members_of_reference=[
                             dict(gene=g, ref_P=arms[ref][g]['P'], test_P=arms[test][g]['P'],
                                  dlog10P=-math.log10(max(arms[test][g]['P'], tiny))
                                          + math.log10(max(arms[ref][g]['P'], tiny)))
                             for g in band[f'{ref}_in_band']])

    C.atomic_tsv(out / f'{args.tag}.per_gene.tsv',
                 ['gene', 'truth_P'] + [f'P_{a}' for a in ARMS] +
                 ['dlog10P_phi_vs_flat', 'dlog10P_phibeta_vs_none', 'MAC', 'Number_rare'],
                 [[g, truth[g]['P']] + [arms[a][g]['P'] for a in ARMS] +
                  [-math.log10(max(arms['phi'][g]['P'], tiny)) + math.log10(max(arms['flat'][g]['P'], tiny)),
                   -math.log10(max(arms['phibeta'][g]['P'], tiny)) + math.log10(max(arms['none'][g]['P'], tiny)),
                   arms['none'][g]['MAC'], arms['none'][g]['nrare']] for g in genes])
    C.atomic_tsv(out / f'{args.tag}.threshold_counts.tsv',
                 ['threshold', 'truth'] + ARMS,
                 [[f'{t:g}', counts['truth'][f'{t:g}']] + [counts[a][f'{t:g}'] for a in ARMS] for t in thr])
    C.atomic_tsv(out / f'{args.tag}.truth_genes.tsv',
                 ['gene', 'truth_P'] + [f'{a}_P' for a in ARMS] +
                 ['dlog10P_phibeta_vs_none', 'phibeta_lower_than_none',
                  'dlog10P_phi_vs_flat', 'phi_lower_than_flat'],
                 [[r['gene'], r['truth_P']] + [r[f'{a}_P'] for a in ARMS] +
                  [r['dlog10P_phibeta_vs_none'], r['phibeta_lower_than_none'],
                   r['dlog10P_phi_vs_flat'], r['phi_lower_than_flat']] for r in trows])
    C.atomic_json(out / f'{args.tag}.summary.json',
                  dict(status='DIRECTION CHECK ONLY -- not a verdict; no permutation null was run',
                       gate=gate, n_genes=len(genes), thresholds=thr,
                       threshold_counts=counts, questions=questions,
                       truth_alpha=args.truth_alpha, n_truth_significant=len(tsig),
                       truth_genes=trows, just_below_band=band,
                       parallelism_used=json.loads(args.par_record),
                       arm_definitions={
                           'none': 'group file var+anno only; SAIGE default Beta(MAF;1,25)',
                           'flat': 'weight line = 1.0 (equals --is_no_weight_in_groupTest=TRUE)',
                           'phi': 'weight line = phi (probability, unchanged); REPLACES the MAF weight',
                           'phibeta': 'weight line = phi * Beta(MAF;1,25); MAF weight retained'}))
    C.atomic_json(out / f'{args.tag}.done', dict(status='complete', n_genes=len(genes)))
    m = questions['mechanism_phi_vs_flat']['dlog10P']
    p = questions['paper_phibeta_vs_none']['dlog10P']
    print('CMP4_DONE genes=%d gate_none_equals_truth=%s' % (len(genes), gate['none_equals_truth']))
    print('mechanism phi>flat frac=%.4f median_dlog10P=%+.4f' % (m['frac_first_more_significant'], m['median']))
    print('paper phibeta>none frac=%.4f median_dlog10P=%+.4f' % (p['frac_first_more_significant'], p['median']))

if __name__ == '__main__':
    main()
