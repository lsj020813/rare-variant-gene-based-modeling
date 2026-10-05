#!/usr/bin/env python3
import argparse
import csv
import gzip
import json
import os
from pathlib import Path

import l3_common as C

def parser():
    p = argparse.ArgumentParser(description='layer-3 materials inventory')
    p.add_argument('--root', required=True, help='band reference root (target band)')
    p.add_argument('--alt-root', required=True, help="comparison band root, or 'none'")
    p.add_argument('--features-dir', required=True)
    p.add_argument('--cadd-extract-dir', required=True)
    p.add_argument('--resid-dir', required=True)
    p.add_argument('--truth-dir', required=True)
    p.add_argument('--gwas-dir', required=True)
    p.add_argument('--gencode', required=True)
    p.add_argument('--band-vcf-dir', required=True)
    p.add_argument('--phi-json', required=True)
    p.add_argument('--traits', required=True)
    p.add_argument('--threads', required=True, type=int)
    p.add_argument('--out', required=True)
    return p

def fileinfo(path):
    p = Path(path)
    if not p.exists():
        return dict(path=str(p), exists=False)
    st = p.stat()
    return dict(path=str(p), exists=True, bytes=st.st_size,
                mtime=int(st.st_mtime), is_dir=p.is_dir())

def main():
    args = parser().parse_args()
    C.load_libraries(args.threads)
    np = C.np
    R = Path(args.root)
    traits = [t.strip() for t in args.traits.split(',') if t.strip()]
    inv = dict(root=str(R))

    zpath = R / 'annot/cache/fm_all.npz'
    z = np.load(zpath, allow_pickle=True)
    X, cols = z['X'], list(z['cols'].astype(str))
    keys, genes, chrs = z['key37'].astype(str), z['gene'].astype(str), z['chr'].astype(int).astype(str)
    maf = X[:, cols.index('maf')].astype(float)
    per_chr = {}
    for ch in sorted(set(chrs), key=int):
        sel = chrs == ch
        per_chr[ch] = dict(pairs=int(sel.sum()), distinct_variants=len(set(keys[sel].tolist())),
                           genes=len(set(genes[sel].tolist())))
    gene_ids = sorted(set(genes.tolist()))
    inv['cache'] = dict(**fileinfo(zpath), shape=list(X.shape), columns=cols,
                        n_pairs=int(X.shape[0]), n_distinct_variants=len(set(keys.tolist())),
                        n_genes=len(gene_ids), chromosomes=sorted(set(chrs), key=int),
                        chrX_present=('X' in set(chrs) or '23' in set(chrs)),
                        maf_min=float(np.nanmin(maf)), maf_max=float(np.nanmax(maf)),
                        gene_id_style=('ENSG' if gene_ids[0].startswith('ENSG') else 'other'),
                        gene_id_has_version=bool('.' in gene_ids[0]),
                        stats_json=fileinfo(R / 'annot/cache/fm_all.stats.json'),
                        na_fraction={c: float(np.isnan(X[:, i]).mean()) for i, c in enumerate(cols)})
    summaries = sorted((R / 'groupfiles_bwg').glob('chr*.summary.json'))
    inv['groupfiles'] = dict(n_files=len(summaries),
                             pairs_total=sum(int(json.load(open(p))['pairs']) for p in summaries),
                             matches_cache_rows=(sum(int(json.load(open(p))['pairs']) for p in summaries)
                                                 == int(X.shape[0])),
                             assignment='3kb TSS union rE2G (B arm)', per_chrom=per_chr)

    if args.alt_root != 'none':
        ap = Path(args.alt_root) / 'annot/cache/fm_all.npz'
        if ap.exists():
            za = np.load(ap, allow_pickle=True)
            Xa, ca = za['X'], list(za['cols'].astype(str))
            ma = Xa[:, ca.index('maf')].astype(float)
            ka = za['key37'].astype(str)
            inv['alt_cache'] = dict(**fileinfo(ap), shape=list(Xa.shape),
                                    maf_min=float(np.nanmin(ma)), maf_max=float(np.nanmax(ma)),
                                    key_overlap_with_target=len(set(ka.tolist()) & set(keys.tolist())))
        else:
            inv['alt_cache'] = fileinfo(ap)

    ds = {}
    for ch in sorted(set(chrs), key=int):
        p = R / f'annot/ds/chr{ch}.ds.npz'
        ds[ch] = fileinfo(p)
    inv['dosage'] = dict(pattern=str(R / 'annot/ds/chr{N}.ds.npz'),
                         n_present=sum(1 for v in ds.values() if v['exists']),
                         total_bytes=sum(v.get('bytes', 0) for v in ds.values()), per_chrom=ds)
    probe = C.Dosage(args.root, sorted(set(chrs), key=int)[-1])
    inv['dosage']['probe'] = dict(chrom=sorted(set(chrs), key=int)[-1], samples=probe.ns,
                                  variants=probe.n_variants, nnz=probe.nnz,
                                  density=probe.nnz / (probe.n_variants * probe.ns))
    del probe

    feat, extr = {}, {}
    for ch in sorted(set(chrs), key=int):
        ck = set(keys[chrs == ch].tolist())
        f = Path(args.features_dir) / f'chr{ch}.features.tsv.gz'
        info = fileinfo(f)
        if info['exists']:
            hit = obs = 0
            vals = []
            with gzip.open(f, 'rt') as fh:
                hdr = fh.readline().rstrip('\n').split('\t')
                j = hdr.index('cadd_phred')
                for line in fh:
                    a = line.rstrip('\n').split('\t')
                    if a[0] in ck:
                        hit += 1
                        if a[j] not in ('', 'NA'):
                            obs += 1
                            vals.append(float(a[j]))
            v = np.array(vals)
            info.update(done_marker=Path(args.features_dir, f'chr{ch}.done').exists(),
                        key_coverage=hit / max(len(ck), 1), value_coverage=obs / max(len(ck), 1),
                        mean=float(v.mean()) if len(v) else None,
                        sd=float(v.std()) if len(v) else None,
                        min=float(v.min()) if len(v) else None, max=float(v.max()) if len(v) else None)
        feat[ch] = info
        e = Path(args.cadd_extract_dir) / f'chr{ch}.cadd.tsv'
        ei = fileinfo(e)
        ei['done_marker'] = Path(str(e) + '.done').exists()
        extr[ch] = ei
    inv['cadd'] = dict(features=dict(pattern=str(Path(args.features_dir) / 'chr{N}.features.tsv.gz'),
                                     column='cadd_phred',
                                     n_complete=sum(1 for v in feat.values() if v.get('done_marker')),
                                     per_chrom=feat),
                       extract=dict(pattern=str(Path(args.cadd_extract_dir) / 'chr{N}.cadd.tsv'),
                                    n_complete=sum(1 for v in extr.values() if v.get('done_marker')),
                                    per_chrom=extr))

    with open(args.phi_json) as fh:
        model = json.load(fh)
    dcols = list(model['design']['cols'])
    inv['phi'] = dict(**fileinfo(args.phi_json), version=model.get('version'),
                      n_design_cols=len(dcols), design_cols=dcols,
                      cols_missing_from_cache=[c for c in dcols if c not in cols],
                      selected=model['design']['selected_columns'],
                      inactive=model['design']['inactive_columns'],
                      order=model['design']['order'], n_coefficients=len(model['coefficients']['shared'])
                      if 'shared' in model['coefficients'] else None,
                      default_intercept=model['default_intercept'],
                      intercept_rule=model['default_intercept_rule'],
                      cadd_train_mu=model['design']['basis']['cadd']['mu'],
                      cadd_train_sd=model['design']['basis']['cadd']['sd'],
                      cadd_train_knot_raw=[model['design']['basis']['cadd']['mu']
                                           + model['design']['basis']['cadd']['sd'] * k
                                           for k in model['design']['basis']['cadd']['knots']],
                      output=model['output'])

    rd = Path(args.resid_dir)
    with open(rd / 'build_resid_v8.summary.json') as fh:
        rs = json.load(fh)
    rinfo = {}
    for t in traits:
        p = rd / f'{t}.resid.tsv'
        i = fileinfo(p)
        if i['exists']:
            with open(p, newline='') as fh:
                rdr = csv.reader(fh, delimiter='\t')
                hdr = next(rdr)
                i['columns'] = hdr
                i['rows'] = sum(1 for _ in rdr)
            m = rs['traits'].get(t, {})
            i.update(manifest_n=m.get('n'), resid_kind=m.get('resid_kind'),
                     covariates=m.get('covariates'), rank=m.get('rank'),
                     sha_matches_manifest=(C.sha_file(p) == m.get('residual_sha256')))
        rinfo[t] = i
    inv['residuals'] = dict(dir=str(rd), builder_version=rs.get('version'),
                            sample_n=rs.get('sample_n'), complete_case_n=rs.get('complete_case_n'),
                            union_n=rs.get('union_n'),
                            done_marker=fileinfo(rd / 'build_resid_v8.done'),
                            residual_policy=rs.get('residual_policy'),
                            n_traits_available=len(list((rd).glob('*.resid.tsv'))), per_trait=rinfo)

    tinfo = {}
    for t in traits:
        present, sig, tested, groups = [], 0, 0, set()
        for ch in sorted(set(chrs), key=int):
            f = Path(args.truth_dir) / f'{t}.chr{ch}'
            if f.is_file() and Path(str(f) + '.done').is_file():
                present.append(ch)
                with open(f, newline='') as fh:
                    for row in csv.DictReader(fh, delimiter='\t'):
                        groups.add((row.get('Group'), row.get('max_MAF')))
                        tested += 1
                        try:
                            if float(row['Pvalue']) < 2.5e-6:
                                sig += 1
                        except (TypeError, ValueError):
                            pass
        tinfo[t] = dict(chromosomes_present=present, n_missing=len(set(chrs)) - len(present),
                        rows=tested, group_maxmaf_values=sorted(str(x) for x in groups),
                        n_sig_at_2p5e6=sig)
    inv['truth'] = dict(dir=args.truth_dir, pattern='<trait>.chr<N> (+ .done)',
                        columns=['Region', 'Group', 'max_MAF', 'Pvalue', 'Pvalue_Burden',
                                 'Pvalue_SKAT', 'BETA_Burden', 'SE_Burden', 'MAC',
                                 'Number_rare', 'Number_ultra_rare'], per_trait=tinfo)
    anyt = next((t for t in traits if tinfo[t]['chromosomes_present']), None)
    if anyt:
        ch = tinfo[anyt]['chromosomes_present'][-1]
        tg = set()
        with open(Path(args.truth_dir) / f'{anyt}.chr{ch}', newline='') as fh:
            for row in csv.DictReader(fh, delimiter='\t'):
                tg.add(row['Region'].split('.')[0])
        cg = set(genes[chrs == ch].tolist())
        inv['truth']['gene_id_join_probe'] = dict(
            chrom=ch, truth_genes=len(tg), cache_genes=len(cg), intersection=len(tg & cg),
            fraction_of_cache_genes=len(tg & cg) / max(len(cg), 1))

    gw = {}
    for t in traits:
        n = sum(1 for ch in sorted(set(chrs), key=int)
                if (Path(args.gwas_dir) / f'{t}.chr{ch}.txt').is_file()
                and (Path(args.gwas_dir) / f'{t}.chr{ch}.txt.done').is_file())
        gw[t] = dict(chromosomes_complete=n, expected=len(set(chrs)))
    hdr = None
    p0 = sorted(Path(args.gwas_dir).glob('*.chr*.txt'))
    if p0:
        with open(p0[0]) as fh:
            hdr = fh.readline().rstrip('\n').split('\t')
    inv['gwas_common'] = dict(dir=args.gwas_dir, pattern='<trait>.chr<N>.txt (+ .done)',
                              columns=hdr, per_trait=gw)
    bv = {ch: fileinfo(Path(args.band_vcf_dir) / f'chr{ch}.band.vcf.gz')
          for ch in sorted(set(chrs), key=int)}
    inv['band_vcf'] = dict(dir=args.band_vcf_dir,
                           n_present=sum(1 for v in bv.values() if v['exists']),
                           total_bytes=sum(v.get('bytes', 0) for v in bv.values()),
                           index_present=all((Path(args.band_vcf_dir) / f'chr{ch}.band.vcf.gz.csi').exists()
                                             for ch in sorted(set(chrs), key=int)),
                           note='dosage is read from the DS npz caches, not the VCFs',
                           per_chrom=bv)
    inv['gencode'] = fileinfo(args.gencode)

    C.atomic_json(Path(args.out) / 'l3_inventory.json', inv)
    C.atomic_json(Path(args.out) / 'l3_inventory.done',
                  dict(status='complete',
                       sha256=C.sha_file(Path(args.out) / 'l3_inventory.json')))
    print('L3_INVENTORY_DONE', flush=True)

if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        msg = str(error) if isinstance(error, RuntimeError) and str(error).startswith('GATE FAIL:') \
            else type(error).__name__
        print(f'[l3_inventory] {msg}', flush=True)
        raise SystemExit(1) from None
