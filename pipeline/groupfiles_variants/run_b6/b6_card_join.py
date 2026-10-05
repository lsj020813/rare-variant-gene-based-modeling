#!/usr/bin/env python3
import os as _config_os
import re as _config_re

def _required(name):
    value = _config_os.environ.get(name)
    if value is None or not value.strip():
        raise RuntimeError(f'Set {name} before running this script')
    return value

def _required_int(name):
    value = int(_required(name))
    if value <= 0:
        raise ValueError(f'{name} must be a positive integer')
    return value

def _configured(value):
    defaults = {'BCFTOOLS': 'bcftools', 'PLINK2': 'plink2', 'PYTHON': 'python3', 'TABIX': 'tabix', 'BGZIP': 'bgzip', 'SAMTOOLS': 'samtools', 'BEDTOOLS': 'bedtools', 'CROSSMAP': 'CrossMap', 'UDOCKER_BIN': 'udocker'}
    def resolve(match):
        name = match.group(1) or match.group(2)
        if name in _config_os.environ:
            return _required(name)
        return defaults.get(name) or _required(name)
    return _config_re.sub(r"\$\{([A-Z][A-Z0-9_]*)\}|\$([A-Z][A-Z0-9_]*)", resolve, value)

def _source_path(base, value):
    return value if _config_os.path.isabs(value) else _config_os.path.join(base, value)
import json, subprocess, sys, os, gzip
import numpy as np
CHROMS = sys.argv[1:] or [f'chr{i}' for i in range(1, 23)]
TMP = _configured('${PROJECT_ROOT}/work/tmp')
LIFT = f'{TMP}/lift38'
CCRE = f'{TMP}/ccre/GRCh38-cCREs.V4.bed'
GNO = f'{TMP}/gnocchi/gnocchi_1kb_qc.txt.gz'
ANN = _configured('${PROJECT_ROOT}/work/ref/annotation_ref/annotation_data')
BCF = _configured('${BCFTOOLS}')
TABIX = _configured('${TABIX}')
BEDT = _configured('${BEDTOOLS}')
OUT = _configured('${PROJECT_ROOT}/work/ref/b6_cards')
os.makedirs(OUT, exist_ok=True)

def nc(c):
    return c[3:] if c.startswith('chr') else c

def sh(cmd, check=True):
    r = subprocess.run(['bash', '-c', cmd], capture_output=True, text=True)
    if check and r.returncode != 0:
        print('FAIL:', cmd[:160], '\n', (r.stderr or r.stdout)[-500:], flush=True)
        sys.exit(1)
    return r
CLASSES = ['PLS', 'pELS', 'dELS', 'CA-CTCF', 'CA-H3K4me3', 'CA-TF', 'CA', 'TF']
if not os.path.exists(f'{OUT}/ccre.s.bed'):
    sh(_configured(f"""awk -F'\\t' '{{print $1"\\t"$2"\\t"$3"\\t"$NF}}' {CCRE} | sort -k1,1 -k2,2n > {OUT}/ccre.s.bed"""))
if not os.path.exists(f'{OUT}/gno.s.bed'):
    sh(f"""zcat {GNO} | tail -n +2 | awk -F'\\t' '{{print $1"\\t"$2"\\t"$3"\\t"$9}}' | sort -k1,1 -k2,2n > {OUT}/gno.s.bed""")
for CN in CHROMS:
    marker = f'{OUT}/{CN}.done'
    if os.path.exists(marker):
        print(f'[{CN}] skip', flush=True)
        continue
    v = f'{LIFT}/{CN}.sites38.vcf.gz'
    assert os.path.exists(v), f'missing {v}'
    rows = []
    p = subprocess.Popen(f"{BCF} query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/MAF\\t%INFO/R2\\n' {v}", shell=True, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for line in p.stdout:
        f = line.rstrip('\n').split('\t')
        try:
            maf, r2 = (float(f[4]), float(f[5]))
        except (ValueError, IndexError):
            continue
        if 0.001 <= maf <= 0.01:
            rows.append((nc(f[0]), int(f[1]), f[2], f[3], maf, r2))
    n = len(rows)
    print(f'[{CN}] band variants: {n:,}', flush=True)
    assert n > 0, f'GATE FAIL {CN}: zero band variants'
    maf = np.array([r[4] for r in rows])
    r2 = np.array([r[5] for r in rows])
    is_indel = np.array([len(r[2]) != 1 or len(r[3]) != 1 for r in rows])
    lo, hi = (rows[0][1], rows[-1][1])
    with open(f'{OUT}/{CN}.var.bed', 'w') as fh:
        for i, r in enumerate(rows):
            fh.write(f'chr{r[0]}\t{r[1] - 1}\t{r[1]}\t{i}\n')
    sh(f'sort -k1,1 -k2,2n {OUT}/{CN}.var.bed > {OUT}/{CN}.var.s.bed')
    Z = np.zeros((n, len(CLASSES)), dtype=np.int8)
    cidx = {c: k for k, c in enumerate(CLASSES)}
    sh(f'{BEDT} intersect -a {OUT}/{CN}.var.s.bed -b {OUT}/ccre.s.bed -wa -wb > {OUT}/{CN}.ccre.tsv || true', check=False)
    if os.path.exists(f'{OUT}/{CN}.ccre.tsv'):
        for line in open(f'{OUT}/{CN}.ccre.tsv'):
            f = line.rstrip('\n').split('\t')
            if len(f) < 8:
                continue
            i = int(f[3])
            for cc in f[7].split(','):
                cc = cc.strip()
                if cc in cidx:
                    Z[i, cidx[cc]] = 1
    nlab = Z.sum(1)
    gz = np.full(n, np.nan)
    sh(f'{BEDT} intersect -a {OUT}/{CN}.var.s.bed -b {OUT}/gno.s.bed -wa -wb > {OUT}/{CN}.gno.tsv || true', check=False)
    if os.path.exists(f'{OUT}/{CN}.gno.tsv'):
        for line in open(f'{OUT}/{CN}.gno.tsv'):
            f = line.rstrip('\n').split('\t')
            if len(f) < 8:
                continue
            try:
                gz[int(f[3])] = float(f[7])
            except ValueError:
                pass
    by = {}
    for i, r in enumerate(rows):
        by.setdefault((r[1], r[2], r[3]), []).append(i)
    REG_BARE = f'{nc(CN)}:{lo}-{hi}'
    cadd = np.full(n, np.nan)
    for f_ in (f'{ANN}/cadd/whole_genome_SNVs.tsv.gz', f'{ANN}/cadd/gnomad.genomes.r3.0.indel.tsv.gz'):
        if not os.path.exists(f_):
            continue
        pr = subprocess.Popen(f'{TABIX} {f_} {REG_BARE} 2>/dev/null', shell=True, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
        for line in pr.stdout:
            f = line.rstrip('\n').split('\t')
            if len(f) < 6:
                continue
            for i in by.get((int(f[1]), f[2], f[3]), ()):
                try:
                    cadd[i] = float(f[5])
                except ValueError:
                    pass
    sp = {k: np.full(n, np.nan) for k in ('ag', 'al', 'dg', 'dl')}
    for f_ in (f'{ANN}/spliceAI/spliceai_scores.raw.snv.hg38.vcf.gz', f'{ANN}/spliceAI/spliceai_scores.raw.indel.hg38.vcf.gz'):
        if not os.path.exists(f_):
            continue
        pr = subprocess.Popen(f"{BCF} query -r {REG_BARE} -f '%POS\\t%REF\\t%ALT\\t%INFO/SpliceAI\\n' {f_} 2>/dev/null", shell=True, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
        for line in pr.stdout:
            f = line.rstrip('\n').split('\t')
            if len(f) < 4:
                continue
            hit = by.get((int(f[0]), f[1], f[2]))
            if not hit:
                continue
            parts = f[3].split('|')
            if len(parts) < 6:
                continue
            try:
                vals = [float(x) for x in parts[2:6]]
            except ValueError:
                continue
            for i in hit:
                for k, val in zip(('ag', 'al', 'dg', 'dl'), vals):
                    sp[k][i] = val
    cards = {'MAF': maf, 'R2': r2, 'is_indel': is_indel.astype(float), 'gnocchi_z': gz, 'cadd_phred': cadd, 'n_ccre_labels': nlab.astype(float), 'spliceai_ag': sp['ag'], 'spliceai_al': sp['al'], 'spliceai_dg': sp['dg'], 'spliceai_dl': sp['dl']}
    for k, cname in enumerate(CLASSES):
        cards[f'lab_{cname}'] = Z[:, k].astype(float)
    gate = {}
    for nm, arr in cards.items():
        fin = np.isfinite(arr)
        cov = float(fin.mean())
        vals = arr[fin]
        uniq = int(len(np.unique(vals))) if vals.size else 0
        if vals.size:
            _, cnts = np.unique(vals, return_counts=True)
            mode_share = float(cnts.max() / vals.size)
        else:
            mode_share = 1.0
        gate[nm] = {'coverage': round(cov, 5), 'n_unique': uniq, 'mode_share': round(mode_share, 5), 'std': round(float(vals.std()), 5) if vals.size else 0.0, 'n_negative': int((vals < 0).sum()) if vals.size else 0, 'min': round(float(vals.min()), 5) if vals.size else None, 'max': round(float(vals.max()), 5) if vals.size else None, 'pass_info': bool(mode_share < 0.9 and uniq > 5 and (cov > 0.3))}
    names = [k for k in cards if gate[k]['coverage'] > 0.05]
    corr, strong = ({}, [])
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = (cards[names[i]], cards[names[j]])
            ok = np.isfinite(a) & np.isfinite(b)
            if ok.sum() < 200:
                continue
            if a[ok].std() == 0 or b[ok].std() == 0:
                continue
            cc = float(np.corrcoef(a[ok], b[ok])[0, 1])
            corr[f'{names[i]}|{names[j]}'] = round(cc, 4)
            if abs(cc) >= 0.7:
                strong.append((names[i], names[j], round(cc, 4)))
    fin = np.isfinite(gz)
    by_class = {}
    for k, cname in enumerate(CLASSES):
        sel = fin & (Z[:, k] == 1)
        if sel.sum() >= 30:
            by_class[cname] = {'n': int(sel.sum()), 'mean_z': round(float(gz[sel].mean()), 4), 'frac_z_gt_2': round(float((gz[sel] > 2).mean()), 4)}
    bg = fin & (nlab == 0)
    background = {'n': int(bg.sum()), 'mean_z': round(float(gz[bg].mean()), 4), 'frac_z_gt_2': round(float((gz[bg] > 2).mean()), 4)} if bg.sum() >= 30 else None
    led = {'chrom': CN, 'n_band': n, 'span': [lo, hi], 'gates': gate, 'corr': corr, 'strong_pairs_ge_0.7': strong, 'ccre_counts': {c: int(Z[:, k].sum()) for k, c in enumerate(CLASSES)}, 'labeled_frac': round(float((nlab > 0).mean()), 5), 'labels_per_variant': {str(k): int((nlab == k).sum()) for k in range(0, min(5, int(nlab.max()) + 1))}, 'gnocchi_by_class': by_class, 'gnocchi_background_unlabeled': background}
    json.dump(led, open(f'{OUT}/{CN}.b6.json', 'w'), indent=1)
    for f_ in (f'{CN}.var.bed', f'{CN}.ccre.tsv', f'{CN}.gno.tsv'):
        try:
            os.remove(f'{OUT}/{f_}')
        except OSError:
            pass
    open(marker, 'w').close()
    print(f'[{CN}] gates: ' + ', '.join((f"{k}={gate[k]['coverage']:.2f}" for k in ('gnocchi_z', 'cadd_phred', 'spliceai_ag', 'n_ccre_labels'))), flush=True)
print('B6_JOIN_COMPLETE')
