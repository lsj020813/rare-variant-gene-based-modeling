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

def load(path):
    d = {}
    with open(path, newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            try:
                d[r['Region']] = dict(P=float(r['Pvalue']),
                                      PB=float(r['Pvalue_Burden']) if r.get('Pvalue_Burden') else None,
                                      PS=float(r['Pvalue_SKAT']) if r.get('Pvalue_SKAT') else None,
                                      MAC=float(r['MAC']) if r.get('MAC') else None,
                                      nrare=r.get('Number_rare'))
            except (TypeError, ValueError):
                continue
    return d

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--none', required=True)
    ap.add_argument('--flat', required=True)
    ap.add_argument('--phi', required=True)
    ap.add_argument('--truth', required=True)
    ap.add_argument('--truth-alpha', required=True, type=float)
    ap.add_argument('--thresholds', required=True)
    ap.add_argument('--just-below', required=True, help='lo,hi for the (d) band, e.g. 2.5e-6,1e-5')
    ap.add_argument('--out', required=True)
    ap.add_argument('--tag', required=True)
    args = ap.parse_args()
    C.load_libraries(1)
    np = C.np
    from scipy.stats import spearmanr, kendalltau
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)

    arms = dict(none=load(args.none), flat=load(args.flat), phi=load(args.phi))
    truth = load(args.truth)
    genes = sorted(set.intersection(*[set(v) for v in arms.values()]))
    C.require(bool(genes), 'no genes shared by the three arms')
    thr = [float(x) for x in args.thresholds.split(',')]
    lo, hi = (float(x) for x in args.just_below.split(','))

    counts = {a: {f'{t:g}': int(sum(arms[a][g]['P'] <= t for g in genes)) for t in thr} for a in arms}

    Pp = np.array([arms['phi'][g]['P'] for g in genes])
    Pf = np.array([arms['flat'][g]['P'] for g in genes])
    Pn = np.array([arms['none'][g]['P'] for g in genes])
    tiny = np.finfo(float).tiny
    d = -np.log10(np.maximum(Pp, tiny)) + np.log10(np.maximum(Pf, tiny))
    sp = spearmanr(Pp, Pf); kt = kendalltau(Pp, Pf)
    rank = dict(spearman_rho=float(sp.statistic), spearman_p=float(sp.pvalue),
                kendall_tau=float(kt.statistic), kendall_p=float(kt.pvalue),
                dlog10P_phi_minus_flat=dict(
                    mean=float(d.mean()), sd=float(d.std()), median=float(np.median(d)),
                    q={f'p{q}': float(np.percentile(d, q)) for q in (1, 5, 25, 50, 75, 95, 99)},
                    frac_phi_more_significant=float((d > 0).mean()),
                    frac_abs_gt_0p5=float((np.abs(d) > 0.5).mean()),
                    frac_abs_gt_1=float((np.abs(d) > 1).mean()),
                    max_gain=float(d.max()), max_loss=float(d.min())),
                spearman_phi_vs_none=float(spearmanr(Pp, Pn).statistic))

    tsig = sorted(g for g in truth if truth[g]['P'] < args.truth_alpha)
    trows = []
    for g in tsig:
        r = dict(gene=g, truth_P=truth[g]['P'])
        for a in ('none', 'flat', 'phi'):
            r[f'{a}_P'] = arms[a][g]['P'] if g in arms[a] else None
        if r['phi_P'] and r['flat_P']:
            r['dlog10P_phi_minus_flat'] = -math.log10(max(r['phi_P'], tiny)) + math.log10(max(r['flat_P'], tiny))
            r['phi_lower_than_flat'] = bool(r['phi_P'] < r['flat_P'])
            r['phi_lower_than_truth'] = bool(r['phi_P'] < r['truth_P'])
        trows.append(r)

    band = dict(lo=lo, hi=hi)
    for a in ('none', 'flat', 'phi'):
        band[f'{a}_in_band'] = sorted(g for g in genes if lo < arms[a][g]['P'] <= hi)
    crossed = sorted(g for g in genes if arms['flat'][g]['P'] > lo >= arms['phi'][g]['P'])
    lost = sorted(g for g in genes if arms['phi'][g]['P'] > lo >= arms['flat'][g]['P'])
    band['crossed_into_significance_under_phi'] = crossed
    band['lost_significance_under_phi'] = lost
    band['n_crossed'] = len(crossed)
    band['n_lost'] = len(lost)

    C.atomic_tsv(out / f'{args.tag}.per_gene.tsv',
                 ['gene', 'P_none', 'P_flat', 'P_phi', 'dlog10P_phi_minus_flat', 'MAC_flat', 'Number_rare_flat'],
                 [[g, arms['none'][g]['P'], arms['flat'][g]['P'], arms['phi'][g]['P'],
                   -math.log10(max(arms['phi'][g]['P'], tiny)) + math.log10(max(arms['flat'][g]['P'], tiny)),
                   arms['flat'][g]['MAC'], arms['flat'][g]['nrare']] for g in genes])
    C.atomic_tsv(out / f'{args.tag}.threshold_counts.tsv',
                 ['threshold'] + list(arms),
                 [[f'{t:g}'] + [counts[a][f'{t:g}'] for a in arms] for t in thr])
    C.atomic_tsv(out / f'{args.tag}.truth_genes.tsv',
                 ['gene', 'truth_P', 'none_P', 'flat_P', 'phi_P', 'dlog10P_phi_minus_flat',
                  'phi_lower_than_flat', 'phi_lower_than_truth'],
                 [[r['gene'], r['truth_P'], r.get('none_P'), r.get('flat_P'), r.get('phi_P'),
                   r.get('dlog10P_phi_minus_flat'), r.get('phi_lower_than_flat'),
                   r.get('phi_lower_than_truth')] for r in trows])
    C.atomic_json(out / f'{args.tag}.summary.json',
                  dict(status='DIRECTION CHECK ONLY -- not a verdict; no permutation null was run',
                       n_genes=len(genes), thresholds=thr, threshold_counts=counts,
                       rank_and_delta=rank, truth_alpha=args.truth_alpha,
                       n_truth_significant=len(tsig), truth_genes=trows, just_below_band=band))
    C.atomic_json(out / f'{args.tag}.done', dict(status='complete', n_genes=len(genes)))
    print('CMP_DONE genes=%d spearman=%.6f frac_phi_better=%.4f n_truth_sig=%d crossed=%d lost=%d'
          % (len(genes), rank['spearman_rho'], rank['dlog10P_phi_minus_flat']['frac_phi_more_significant'],
             len(tsig), len(crossed), len(lost)))

if __name__ == '__main__':
    main()
