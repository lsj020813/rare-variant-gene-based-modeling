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

ARMS = ['none', 'flat', 'phi', 'phibeta', 'phisqrt', 'permA1', 'permA2', 'permA3']
PERMA = ['permA1', 'permA2', 'permA3']
PCOLS = [('P', 'Pvalue'), ('PB', 'Pvalue_Burden'), ('PS', 'Pvalue_SKAT')]

def load(path):
    d = {}
    for r in csv.DictReader(open(path), delimiter='\t'):
        rec = {}
        for short, col in PCOLS:
            try:
                rec[short] = float(r[col])
            except (TypeError, ValueError, KeyError):
                rec[short] = None
        rec['MAC'], rec['nrare'] = r.get('MAC'), r.get('Number_rare')
        if rec['P'] is not None:
            d[r['Region']] = rec
    return d

def signtest(np, A, B, genes, key):
    a = np.array([A[g][key] for g in genes if A[g][key] is not None and B[g][key] is not None])
    b = np.array([B[g][key] for g in genes if A[g][key] is not None and B[g][key] is not None])
    tiny = np.finfo(float).tiny
    d = -np.log10(np.maximum(a, tiny)) + np.log10(np.maximum(b, tiny))
    nz = d[d != 0]
    k, n = int((nz > 0).sum()), int(len(nz))
    z = (k - n / 2) / math.sqrt(n / 4) if n else None
    return dict(n_compared=int(len(d)), n_nontied=n, n_test_better=k,
                frac_test_better=(k / n if n else None),
                median_dlog10P=float(np.median(d)), mean_dlog10P=float(d.mean()),
                sign_z=(float(z) if z is not None else None),
                max_gain=float(d.max()), max_loss=float(d.min()))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', required=True, help='directory holding <arm>.merged')
    ap.add_argument('--truth', required=True)
    ap.add_argument('--truth-alpha', required=True, type=float)
    ap.add_argument('--thresholds', required=True)
    ap.add_argument('--track-genes', required=True)
    ap.add_argument('--par-record', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--tag', required=True)
    a = ap.parse_args()
    C.load_libraries(1)
    np = C.np
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    tiny = np.finfo(float).tiny

    arms = {x: load(Path(a.dir) / f'{x}.merged') for x in ARMS}
    truth = load(a.truth)
    genes = sorted(set.intersection(*[set(v) for v in arms.values()]) & set(truth))
    C.require(bool(genes), 'no genes shared by all arms and the truth file')
    thr = [float(t) for t in a.thresholds.split(',')]
    extra = [g.strip() for g in a.track_genes.split(',') if g.strip()]

    gate = dict(genes_compared=len(genes),
                none_equals_truth=all(arms['none'][g]['P'] == truth[g]['P'] for g in genes),
                worst_abs_diff=max(abs(arms['none'][g]['P'] - truth[g]['P']) for g in genes))

    counts = {}
    for key, _ in PCOLS:
        counts[key] = {x: {f'{t:g}': int(sum(arms[x][g][key] is not None and arms[x][g][key] <= t
                                             for g in genes)) for t in thr} for x in ARMS}
        counts[key]['truth'] = {f'{t:g}': int(sum(truth[g][key] is not None and truth[g][key] <= t
                                                  for g in genes)) for t in thr}
        cs = [counts[key][p] for p in PERMA]
        counts[key]['permA_mean'] = {f'{t:g}': float(np.mean([c[f'{t:g}'] for c in cs])) for t in thr}
        counts[key]['permA_min'] = {f'{t:g}': int(min(c[f'{t:g}'] for c in cs)) for t in thr}
        counts[key]['permA_max'] = {f'{t:g}': int(max(c[f'{t:g}'] for c in cs)) for t in thr}

    contrasts = {}
    for x in ARMS:
        for ref in ('flat', 'none'):
            if x == ref:
                continue
            for key, col in PCOLS[:1] + PCOLS[2:]:
                contrasts[f'{x}|vs_{ref}|{col}'] = signtest(np, arms[x], arms[ref], genes, key)
    for key, col in PCOLS[:1] + PCOLS[2:]:
        contrasts[f'phisqrt|vs_permA_pooled|{col}'] = dict(
            note='phisqrt sign_z minus the permA sign_z, same reference',
            vs_flat_phisqrt=contrasts[f'phisqrt|vs_flat|{col}']['sign_z'],
            vs_flat_permA=[contrasts[f'{p}|vs_flat|{col}']['sign_z'] for p in PERMA],
            vs_none_phisqrt=contrasts[f'phisqrt|vs_none|{col}']['sign_z'],
            vs_none_permA=[contrasts[f'{p}|vs_none|{col}']['sign_z'] for p in PERMA])

    tsig = sorted((g for g in genes if truth[g]['P'] < a.truth_alpha), key=lambda g: truth[g]['P'])
    tracked = tsig + [g for g in extra if g in genes and g not in tsig]
    thdr = ['gene', 'group', 'truth_P'] + [f'{x}_P' for x in ARMS] + \
           ['permA_P_mean', 'permA_P_min', 'permA_P_max', 'truth_P_SKAT', 'phisqrt_P_SKAT']
    trows = []
    for g in tracked:
        pv = [arms[p][g]['P'] for p in PERMA]
        trows.append([g, ('truth_significant' if g in tsig else 'just_below'), truth[g]['P']] +
                     [arms[x][g]['P'] for x in ARMS] +
                     [float(np.mean(pv)), float(min(pv)), float(max(pv)),
                      truth[g]['PS'], arms['phisqrt'][g]['PS']])
    C.atomic_tsv(out / f'{a.tag}.tracked_genes.tsv', thdr, trows)

    crows = []
    for key, col in PCOLS:
        for t in thr:
            crows.append([col, f'{t:g}', counts[key]['truth'][f'{t:g}']] +
                         [counts[key][x][f'{t:g}'] for x in ARMS] +
                         [counts[key]['permA_mean'][f'{t:g}'], counts[key]['permA_min'][f'{t:g}'],
                          counts[key]['permA_max'][f'{t:g}']])
    C.atomic_tsv(out / f'{a.tag}.threshold_counts.tsv',
                 ['pvalue_column', 'threshold', 'truth'] + ARMS +
                 ['permA_mean', 'permA_min', 'permA_max'], crows)

    srows = []
    for k, v in sorted(contrasts.items()):
        if 'pooled' in k:
            continue
        arm, ref, col = k.split('|')
        srows.append([arm, ref.replace('vs_', ''), col, v['n_nontied'], v['median_dlog10P'],
                      v['frac_test_better'], v['sign_z'], v['max_gain'], v['max_loss']])
    C.atomic_tsv(out / f'{a.tag}.sign_tests.tsv',
                 ['arm', 'reference', 'pvalue_column', 'n_nontied', 'median_dlog10P',
                  'frac_test_better', 'sign_z', 'max_gain', 'max_loss'], srows)

    hdr = ['gene', 'truth_P'] + [f'{x}_P' for x in ARMS] + \
          [f'{x}_P_SKAT' for x in ARMS] + ['MAC', 'Number_rare']
    C.atomic_tsv(out / f'{a.tag}.per_gene.tsv', hdr,
                 [[g, truth[g]['P']] + [arms[x][g]['P'] for x in ARMS] +
                  [arms[x][g]['PS'] for x in ARMS] +
                  [arms['none'][g]['MAC'], arms['none'][g]['nrare']] for g in genes])

    C.atomic_json(out / f'{a.tag}.summary.json',
                  dict(status='DIRECTION CHECK ONLY -- not a verdict',
                       gate=gate, n_genes=len(genes), thresholds=thr,
                       threshold_counts=counts, contrasts=contrasts,
                       truth_alpha=a.truth_alpha, n_truth_significant=len(tsig),
                       tracked_genes=[dict(zip(thdr, r)) for r in trows],
                       parallelism_used=json.loads(a.par_record),
                       sign_test_definition='z = (k - n/2)/sqrt(n/4) over genes with a non-zero '
                                            'difference in -log10 P; k = test arm more significant',
                       permutation_null_kind='annotation permutation (within gene), NOT a '
                                             'phenotype permutation'))
    C.atomic_json(out / f'{a.tag}.done', dict(status='complete', n_genes=len(genes)))

    print('FINAL8_DONE genes=%d gate=%s' % (len(genes), gate['none_equals_truth']))
    for x in ARMS:
        o = counts['P'][x]['2.5e-06'] if '2.5e-06' in counts['P'][x] else counts['P'][x][f'{thr[0]:g}']
        s = counts['PS'][x][f'{thr[0]:g}']
        zf = contrasts.get(f'{x}|vs_flat|Pvalue', {}).get('sign_z')
        zn = contrasts.get(f'{x}|vs_none|Pvalue', {}).get('sign_z')
        print('%-8s omni=%2d skat=%2d  z_vs_flat=%s  z_vs_none=%s'
              % (x, o, s, ('%+.2f' % zf) if zf is not None else '   . ',
                 ('%+.2f' % zn) if zn is not None else '   . '))

if __name__ == '__main__':
    main()
