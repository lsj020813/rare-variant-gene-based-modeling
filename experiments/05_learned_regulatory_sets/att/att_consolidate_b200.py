import json
import hashlib
import numpy as np
import att_common_b200 as C

def main():
    C.ensure_dirs()
    d = C.master()
    dom_chrom = np.asarray(d['dom_chrom'])
    chs = sorted(set(x.replace('chr', '') for x in dom_chrom), key=int)
    from numpy.lib.format import open_memmap
    cols = [open_memmap(C.BURD + '/col%02d.npy' % i, mode='w+',
                        dtype=np.float32, shape=(C.NDOM, C.NTOT))
            for i in range(C.NCOL)]
    seen = np.zeros(C.NDOM, dtype=bool)
    rep = []
    for ch in chs:
        a = np.load(C.BURD + '/burden_chr%s.npy' % ch, mmap_mode='r')
        gs = np.load(C.BURD + '/domidx_chr%s.npy' % ch)
        for ii, g in enumerate(gs):
            blk = np.asarray(a[ii])
            for cidx in range(C.NCOL):
                cols[cidx][g, :] = blk[cidx]
            seen[g] = True
        rep += json.load(open(C.ATT + '/burden_chr%s.json' % ch))['per_domain']
        del a
    assert seen.all(), '누락 도메인 %d' % int((~seen).sum())
    h0 = hashlib.sha256()
    for cidx in (C.C_TOT_ALL, C.C_CCRE, C.C_MAJ0):
        h0.update(np.ascontiguousarray(cols[cidx][:, :50]).tobytes())
    for cc_ in cols:
        cc_.flush()
    del cols
    Kall = np.load(C.BURD + '/col%02d.npy' % C.C_K_ALL, mmap_mode='r')
    Klab = np.load(C.BURD + '/col%02d.npy' % C.C_K_LAB, mmap_mode='r')
    Kmaj = np.load(C.BURD + '/col%02d.npy' % C.C_K_MAJ, mmap_mode='r')
    Kcc = np.load(C.BURD + '/col%02d.npy' % C.C_K_CCRE, mmap_mode='r')
    f0maj = np.zeros(C.NDOM)
    f0min = np.zeros(C.NDOM)
    fboth = np.zeros(C.NDOM)
    f0all = np.zeros(C.NDOM)
    f0cc = np.zeros(C.NDOM)
    kmean_maj = np.zeros(C.NDOM)
    kmean_min = np.zeros(C.NDOM)
    for g in range(C.NDOM):
        a = np.asarray(Kmaj[g])
        bb = np.asarray(Klab[g]) - a
        f0maj[g] = float((a < 0.5).mean())
        f0min[g] = float((bb < 0.5).mean())
        fboth[g] = float(((a >= 0.5) & (bb >= 0.5)).mean())
        f0all[g] = float((np.asarray(Kall[g]) < 0.5).mean())
        f0cc[g] = float((np.asarray(Kcc[g]) < 0.5).mean())
        kmean_maj[g] = float(a.mean())
        kmean_min[g] = float(bb.mean())
    with open(C.ATT + '/module_carrier_full.csv', 'w') as fh:
        fh.write('domain_idx,frac_K0_major,frac_K0_minor,'
                 'frac_both_modules_K_ge1,frac_K0_domain_all,'
                 'frac_K0_ccre,mean_K_major,mean_K_minor\n')
        for g in range(C.NDOM):
            fh.write('%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.4f,%.4f\n'
                     % (g, f0maj[g], f0min[g], fboth[g], f0all[g], f0cc[g],
                        kmean_maj[g], kmean_min[g]))
    carrier = dict(
        definition='K = DS>=0.5 인 모듈 내 변이 개수, all configured participants',
        n_individuals=int(C.NTOT),
        frac_K0_major_median=float(np.median(f0maj)),
        frac_K0_major_q1=float(np.quantile(f0maj, .25)),
        frac_K0_major_q3=float(np.quantile(f0maj, .75)),
        frac_K0_minor_median=float(np.median(f0min)),
        frac_K0_minor_q1=float(np.quantile(f0min, .25)),
        frac_K0_minor_q3=float(np.quantile(f0min, .75)),
        n_modules_major_all_carriers=int((f0maj == 0).sum()),
        n_modules_minor_all_carriers=int((f0min == 0).sum()),
        frac_modules_all_carriers=float(
            ((f0maj == 0).sum() + (f0min == 0).sum()) / (2.0 * C.NDOM)),
        person_module_pair_frac_K_ge1=float(
            1.0 - 0.5 * (f0maj.mean() + f0min.mean())),
        frac_both_modules_K_ge1_median=float(np.median(fboth)),
        frac_both_modules_K_ge1_mean=float(fboth.mean()),
        frac_K0_domain_all_median=float(np.median(f0all)),
        frac_K0_ccre_median=float(np.median(f0cc)),
        mean_K_major_median=float(np.median(kmean_maj)),
        mean_K_minor_median=float(np.median(kmean_min)))
    del Kall, Klab, Kmaj, Kcc
    mr = np.array([r['match_rate'] for r in rep])
    nm = np.array([r['n_matched'] for r in rep])
    man = dict(
        stage='module_burden', post_hoc=True, amendment='A5',
        n_domains=int(C.NDOM), n_individuals=int(C.NTOT),
        n_variants_universe=int(len(d['key'])),
        maxMAF_main=C.MAF_MAIN, maf_secondary_stratum=C.MAF_SEC,
        n_label_sets=C.NSHUF + 1, seed=C.SEED,
        match_rate_median=float(np.median(mr)),
        match_rate_min=float(mr.min()),
        n_domains_match_rate_lt_099=int((mr < 0.99).sum()),
        n_matched_total=int(nm.sum()),
        n_domains_status_ok=int(sum(1 for r in rep if r['status'] == 'ok')),
        n_domains_ccre_zero=int(sum(1 for r in rep if r['n_ccre'] == 0)),
        burden_checksum_sha256_head=h0.hexdigest(),
        carrier_full_sample=carrier,
        columns={k: v for k, v in dict(
            total_all=C.C_TOT_ALL, total_labelled=C.C_TOT_LAB,
            ccre=C.C_CCRE, total001_all=C.C_TOT001_ALL,
            total001_labelled=C.C_TOT001_LAB, ccre001=C.C_CCRE001,
            maj001=C.C_MAJ001, major_labelset0=C.C_MAJ0,
            K_all=C.C_K_ALL, K_labelled=C.C_K_LAB, K_major=C.C_K_MAJ,
            K_ccre=C.C_K_CCRE).items()},
        note=('개인 수준 값은 fset/att/private/ 에만 존재하며 반출하지 않는다. '
              '열 인덱스는 VCF 순서이고 표본 ID 는 생성하지 않았다.'))
    with open(C.ATT + '/module_burden_manifest.json', 'w') as f:
        json.dump(man, f, indent=1, ensure_ascii=False)
    with open(C.ATT + '/consolidate.done', 'w') as f:
        f.write('ok\n')
    print(json.dumps({k: v for k, v in man.items() if k != 'columns'},
                     indent=1, ensure_ascii=False)[:3000])

if __name__ == '__main__':
    main()
