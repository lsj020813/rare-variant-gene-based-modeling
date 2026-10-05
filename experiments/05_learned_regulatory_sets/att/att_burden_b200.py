import os
import sys
import json
import time
import numpy as np
from cyvcf2 import VCF

import att_common_b200 as C

def main():
    ch = sys.argv[1]
    C.ensure_dirs()
    d = C.master()
    di = np.asarray(d['dom_idx'])
    key = np.asarray(d['key'])
    pos19 = np.asarray(d['pos19'])
    maf = np.asarray(d['maf'], dtype=np.float64)
    hcc = np.asarray(d['has_ccre'], dtype=np.float64)
    dom_chrom = np.asarray(d['dom_chrom'])
    Z = np.load(C.ATT + '/att_labels.npz', allow_pickle=True)
    LAB = Z['labels']
    ok = Z['dom_ok']
    doms = [g for g in range(C.NDOM)
            if dom_chrom[g].replace('chr', '') == str(ch)]
    vcf = VCF(C.VCFD + '/chr%s.vcf.gz' % ch)
    assert len(vcf.samples) == C.NTOT, 'sample 수 불일치 %d' % len(vcf.samples)
    from numpy.lib.format import open_memmap
    out = open_memmap(C.BURD + '/burden_chr%s.npy' % ch, mode='w+',
                      dtype=np.float32, shape=(len(doms), C.NCOL, C.NTOT))
    rep = []
    t0 = time.time()
    for ii, g in enumerate(doms):
        rid = np.where(di == g)[0]
        M = len(rid)
        want = {}
        for j, r in enumerate(rid):
            want.setdefault(key[r], j)
        D = np.zeros((C.NTOT, M), dtype=np.float32)
        got = np.zeros(M, dtype=bool)
        sp = np.sort(np.unique(pos19[rid]))
        iv = []
        a0 = b0 = int(sp[0])
        for x in sp[1:]:
            x = int(x)
            if x - b0 <= 100:
                b0 = x
            else:
                iv.append((a0, b0))
                a0 = b0 = x
        iv.append((a0, b0))
        recs = 0
        for (qa, qb) in iv:
            for v in vcf('%s:%d-%d' % (ch, qa, qb)):
                recs += 1
                for a in v.ALT:
                    k = 'chr%s:%d:%s:%s' % (ch, v.POS, v.REF, a)
                    j = want.get(k)
                    if j is None or got[j]:
                        continue
                    try:
                        ds = np.asarray(v.format('DS'),
                                        dtype=np.float32).ravel()
                    except Exception:
                        ds = None
                    if ds is None or ds.shape[0] != C.NTOT:
                        continue
                    bad = ~np.isfinite(ds) | (ds < 0)
                    if bad.all():
                        continue
                    if bad.any():
                        ds = np.where(bad, np.nanmean(ds[~bad]), ds)
                    D[:, j] = ds
                    got[j] = True
        mrate = float(got.mean())
        W = np.zeros((M, C.NCOL), dtype=np.float32)
        mm = maf[rid]
        cc = hcc[rid] > 0
        lo = mm < C.MAF_SEC
        lab0 = LAB[0][rid]
        W[:, C.C_TOT_ALL] = 1.0
        W[:, C.C_TOT_LAB] = (lab0 >= 0)
        W[:, C.C_CCRE] = cc
        W[:, C.C_TOT001_ALL] = lo
        W[:, C.C_TOT001_LAB] = lo & (lab0 >= 0)
        W[:, C.C_CCRE001] = lo & cc
        W[:, C.C_MAJ001] = lo & (lab0 == 0)
        for L in range(C.NSHUF + 1):
            W[:, C.C_MAJ0 + L] = (LAB[L][rid] == 0)
        W[~got, :] = 0.0
        out[ii, :C.C_K_ALL] = (D @ W[:, :C.C_K_ALL]).T
        Kw = np.zeros((M, 4), dtype=np.float32)
        Kw[:, 0] = 1.0
        Kw[:, 1] = (lab0 >= 0)
        Kw[:, 2] = (lab0 == 0)
        Kw[:, 3] = cc
        Kw[~got, :] = 0.0
        out[ii, C.C_K_ALL:] = ((D >= C.CARRIER_DS).astype(np.float32)
                               @ Kw).T
        del D
        rep.append(dict(domain_idx=int(g), chrom=str(ch), n_var=int(M),
                        n_matched=int(got.sum()), match_rate=mrate,
                        n_records_parsed=int(recs), n_intervals=len(iv),
                        dom_ok=bool(ok[g]),
                        n_ccre=int((cc & got).sum()),
                        n_maf_lt001=int((lo & got).sum()),
                        status='ok' if got.sum() >= 2 else 'too_few_matched'))
    out.flush()
    del out
    np.save(C.BURD + '/domidx_chr%s.npy' % ch,
            np.array(doms, dtype=np.int32))
    with open(C.ATT + '/burden_chr%s.json' % ch, 'w') as f:
        json.dump(dict(chrom=str(ch), n_domains=len(doms),
                       wall_s=round(time.time() - t0, 1), per_domain=rep), f)
    with open(C.ATT + '/burden_chr%s.done' % ch, 'w') as f:
        f.write('%d domains %.1f s\n' % (len(doms), time.time() - t0))

if __name__ == '__main__':
    main()
