import os
import sys
import json
import time
import numpy as np
from cyvcf2 import VCF

import att_common_std as C

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
    Z = np.load(C.LABELS, allow_pickle=True)
    LAB = Z['labels']
    ok = Z['dom_ok']
    assert LAB.shape[0] == C.NPART
    doms = [g for g in range(C.NDOM)
            if dom_chrom[g].replace('chr', '') == str(ch)
            and g % C.DOMSTEP == 0]
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
        sd = D.std(0)
        usable = got & (sd > 1e-6)
        inv = np.zeros(M, dtype=np.float32)
        inv[usable] = 1.0 / sd[usable]
        mm = maf[rid]
        cc = hcc[rid] > 0
        lo = mm < C.MAF_LOW
        W = np.zeros((M, C.NCOL), dtype=np.float32)
        W[:, C.C_RAW_ALL] = got.astype(np.float32)
        W[:, C.C_STD_ALL] = inv
        W[:, C.C_STD_COMMON] = inv * (mm >= C.MAF_COMMON)
        W[:, C.C_STD_LOW] = inv * ((mm >= C.MAF_LOW) & (mm < C.MAF_COMMON))
        W[:, C.C_STD_RARE] = inv * ((mm >= C.MAF_RARE) & (mm < C.MAF_LOW))
        W[:, C.C_STD_VRARE] = inv * (mm < C.MAF_RARE)
        W[:, C.C_RAW_COMMON] = got * (mm >= C.MAF_COMMON)
        W[:, C.C_STD_CCRE] = inv * cc
        W[:, C.C_STD_NOCCRE] = inv * (~cc)
        for L in range(C.NPART):
            lp = LAB[L][rid]
            W[:, C.C_P0 + 2 * L] = inv * (lp == 0)
            W[:, C.C_P0 + 2 * L + 1] = inv * (lp == 1)
        W[~got, :] = 0.0
        out[ii] = (D @ W).T
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
