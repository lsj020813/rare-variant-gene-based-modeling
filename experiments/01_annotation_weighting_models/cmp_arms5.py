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

ARMS = ['none', 'flat', 'phi', 'phibeta', 'phisqrt']
PCOLS = [('P', 'Pvalue'), ('PB', 'Pvalue_Burden'), ('PS', 'Pvalue_SKAT')]

def load(path):
    d = {}
    with open(path, newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            rec = {}
            for short, col in PCOLS:
                v = r.get(col)
                try:
                    rec[short] = float(v)
                except (TypeError, ValueError):
                    rec[short] = None
            rec['MAC'] = r.get('MAC')
            rec['nrare'] = r.get('Number_rare')
            if rec['P'] is not None:
                d[r['Region']] = rec
    return d

def dlog(a, b, tiny):
    if a is None or b is None:
        return None
    return -math.log10(max(a, tiny)) + math.log10(max(b, tiny))

def pair_stats(np, spearmanr, kendalltau, A, B, genes, tiny, key):
    a = np.array([A[g][key] for g in genes if A[g][key] is not None and B[g][key] is not None])
    b = np.array([B[g][key] for g in genes if A[g][key] is not None and B[g][key] is not None])
    if len(a) < 3:
        return dict(n=int(len(a)), note='too few comparable genes')
    d = -np.log10(np.maximum(a, tiny)) + np.log10(np.maximum(b, tiny))
    sp, kt = spearmanr(a, b), kendalltau(a, b)
    return dict(n=int(len(a)), spearman_rho=float(sp.statistic), spearman_p=float(sp.pvalue),
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
    ap.add_argument('--just-below', required=True, help='lo,hi -- the just-below band')
    ap.add_argument('--track-genes', required=True, help='comma list of extra genes to track individually')
    ap.add_argument('--par-record', required=True, help='JSON of arm -> xargs parallelism used')
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
    C.require(bool(genes), 'no genes shared by the five arms and the truth file')
    thr = [float(x) for x in args.thresholds.split(',')]
    lo, hi = (float(x) for x in args.just_below.split(','))
    extra = [g.strip() for g in args.track_genes.split(',') if g.strip()]

    gate = dict(genes_compared=len(genes),
                none_equals_truth=all(arms['none'][g]['P'] == truth[g]['P'] for g in genes),
                worst_abs_diff=max(abs(arms['none'][g]['P'] - truth[g]['P']) for g in genes))

    counts = {}
    for key, _ in PCOLS:
        counts[key] = {a: {f'{t:g}': int(sum(arms[a][g][key] is not None and arms[a][g][key] <= t
                                             for g in genes)) for t in thr} for a in ARMS}
        counts[key]['truth'] = {f'{t:g}': int(sum(truth[g][key] is not None and truth[g][key] <= t
                                                  for g in genes)) for t in thr}

    pairs = [('mechanism_phi_vs_flat', 'phi', 'flat'),
             ('paper_phibeta_vs_none', 'phibeta', 'none'),
             ('sqrt_phisqrt_vs_none', 'phisqrt', 'none'),
             ('scale_phisqrt_vs_phibeta', 'phisqrt', 'phibeta'),
             ('aux_none_vs_flat', 'none', 'flat')]
    questions = {}
    for name, A, B in pairs:
        questions[name] = {k: pair_stats(np, spearmanr, kendalltau, arms[A], arms[B], genes, tiny, k)
                           for k, _ in PCOLS}

    tsig = sorted((g for g in genes if truth[g]['P'] < args.truth_alpha), key=lambda g: truth[g]['P'])
    tracked = tsig + [g for g in extra if g in genes and g not in tsig]
    missing_extra = [g for g in extra if g not in genes]
    trows = []
    for g in tracked:
        r = dict(gene=g, group=('truth_significant' if g in tsig else 'just_below'),
                 truth_P=truth[g]['P'], truth_P_Burden=truth[g]['PB'], truth_P_SKAT=truth[g]['PS'])
        for a in ARMS:
            for k, _ in PCOLS:
                r[f'{a}_{k}'] = arms[a][g][k]
        r['dlog10P_phibeta_vs_none'] = dlog(r['phibeta_P'], r['none_P'], tiny)
        r['dlog10P_phisqrt_vs_none'] = dlog(r['phisqrt_P'], r['none_P'], tiny)
        r['dlog10P_phi_vs_flat'] = dlog(r['phi_P'], r['flat_P'], tiny)
        trows.append(r)

    band = dict(lo=lo, hi=hi, definition=f'{lo:g} < P <= {hi:g} in the reference arm')
    jb_ref = sorted((g for g in genes if lo < truth[g]['P'] <= hi), key=lambda g: truth[g]['P'])
    band['truth_in_band'] = jb_ref
    band['tracked_just_below'] = [g for g in extra if g in genes]
    for name, A, B in pairs[:4]:
        test, ref = A, B
        crossed = sorted(g for g in genes
                         if arms[ref][g]['P'] > lo >= arms[test][g]['P'])
        lostg = sorted(g for g in genes
                       if arms[test][g]['P'] > lo >= arms[ref][g]['P'])
        band[name] = dict(crossed_into_significance=crossed, n_crossed=len(crossed),
                          lost_significance=lostg, n_lost=len(lostg))
    jb = [g for g in extra if g in genes]
    reading = {}
    for a in ARMS:
        n = sum(1 for g in jb if arms[a][g]['P'] is not None and arms[a][g]['P'] < lo)
        reading[a] = dict(n_below_lo=n, of=len(jb),
                          label=('S' if n >= 2 else ('A' if n == 1 else 'direction_only')),
                          per_gene={g: arms[a][g]['P'] for g in jb})
    band['prereg_reading_on_tracked_just_below'] = reading
    band['reading_rule'] = (f'of the {len(jb)} tracked just-below genes: >=2 with P < {lo:g} -> S, '
                            f'1 -> A, 0 -> direction only. No permutation null was run, so every '
                            f'label here is a direction check, not a verdict.')

    hdr = ['gene', 'truth_P', 'truth_P_Burden', 'truth_P_SKAT']
    for a in ARMS:
        hdr += [f'{a}_P', f'{a}_P_Burden', f'{a}_P_SKAT']
    hdr += ['dlog10P_phi_vs_flat', 'dlog10P_phibeta_vs_none', 'dlog10P_phisqrt_vs_none',
            'MAC', 'Number_rare']
    rows = []
    for g in genes:
        r = [g, truth[g]['P'], truth[g]['PB'], truth[g]['PS']]
        for a in ARMS:
            r += [arms[a][g]['P'], arms[a][g]['PB'], arms[a][g]['PS']]
        r += [dlog(arms['phi'][g]['P'], arms['flat'][g]['P'], tiny),
              dlog(arms['phibeta'][g]['P'], arms['none'][g]['P'], tiny),
              dlog(arms['phisqrt'][g]['P'], arms['none'][g]['P'], tiny),
              arms['none'][g]['MAC'], arms['none'][g]['nrare']]
        rows.append(r)
    C.atomic_tsv(out / f'{args.tag}.per_gene.tsv', hdr, rows)

    crows = []
    for key, col in PCOLS:
        for t in thr:
            crows.append([col, f'{t:g}', counts[key]['truth'][f'{t:g}']] +
                         [counts[key][a][f'{t:g}'] for a in ARMS])
    C.atomic_tsv(out / f'{args.tag}.threshold_counts.tsv',
                 ['pvalue_column', 'threshold', 'truth'] + ARMS, crows)

    thdr = ['gene', 'group', 'truth_P', 'truth_P_Burden', 'truth_P_SKAT']
    for a in ARMS:
        thdr += [f'{a}_P', f'{a}_PB', f'{a}_PS']
    thdr += ['dlog10P_phibeta_vs_none', 'dlog10P_phisqrt_vs_none', 'dlog10P_phi_vs_flat']
    C.atomic_tsv(out / f'{args.tag}.tracked_genes.tsv', thdr,
                 [[r[h] for h in thdr] for r in trows])

    C.atomic_json(out / f'{args.tag}.summary.json',
                  dict(status='DIRECTION CHECK ONLY -- not a verdict; no permutation null was run',
                       gate=gate, n_genes=len(genes), thresholds=thr,
                       threshold_counts=counts, questions=questions,
                       truth_alpha=args.truth_alpha, n_truth_significant=len(tsig),
                       tracked_genes=trows, tracked_extra_not_found=missing_extra,
                       just_below_band=band,
                       parallelism_used=json.loads(args.par_record),
                       arm_definitions={
                           'none': 'var+anno only; SAIGE default Beta(MAF;1,25)',
                           'flat': 'weight 1.0 (equals --is_no_weight_in_groupTest=TRUE)',
                           'phi': 'weight phi (probability, unchanged); REPLACES the MAF weight',
                           'phibeta': 'weight phi * Beta(MAF;1,25)',
                           'phisqrt': 'weight sqrt(phi) * Beta(MAF;1,25); exponent fixed at 0.5'}))
    C.atomic_json(out / f'{args.tag}.done', dict(status='complete', n_genes=len(genes)))

    print('CMP5_DONE genes=%d gate_none_equals_truth=%s' % (len(genes), gate['none_equals_truth']))
    for name, _, _ in pairs[:4]:
        q = questions[name]['P']
        s = questions[name]['PS']
        print('%-26s P: frac_first_better=%.4f med_dlog=%+.4f | SKAT: frac=%.4f med_dlog=%+.4f'
              % (name, q['dlog10P']['frac_first_more_significant'], q['dlog10P']['median'],
                 s['dlog10P']['frac_first_more_significant'], s['dlog10P']['median']))
    for a in ARMS:
        print('reading %-8s %s (%d/%d below %g)' % (a, reading[a]['label'], reading[a]['n_below_lo'],
                                                    reading[a]['of'], lo))

if __name__ == '__main__':
    main()
