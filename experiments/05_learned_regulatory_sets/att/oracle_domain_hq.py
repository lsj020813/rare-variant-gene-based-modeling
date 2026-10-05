import os
import sys
import json
import time
import numpy as np
from cyvcf2 import VCF

import att_common_oracle_hq as C

def main():
    ch = sys.argv[1]
    C.ensure_dirs()
    d = C.master()
    di = np.asarray(d['dom_idx'])
    key = np.asarray(d['key'])
    pos19 = np.asarray(d['pos19'])
    maf = np.asarray(d['maf'], dtype=np.float64)
    r2 = np.asarray(d['r2'], dtype=np.float64)
    hcc = np.asarray(d['has_ccre'], dtype=np.float64)
    dom_chrom = np.asarray(d['dom_chrom'])
    Z = np.load(C.LABELS, allow_pickle=True)
    ok = Z['dom_ok']
    PH = C.pheno()
    doms = [g for g in range(C.NDOM)
            if dom_chrom[g].replace('chr', '') == str(ch)
            and g % C.DOMSTEP == 0]
    vcf = VCF(C.VCFD + '/chr%s.vcf.gz' % ch)
    assert len(vcf.samples) == C.NTOT, 'sample 수 불일치 %d' % len(vcf.samples)
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
        mm = maf[rid]
        rep.append(C.oracle_domain(D, got, mm, r2[rid], int(g), str(ch), PH))
        rep[-1].update(dict(n_matched=int(got.sum()), match_rate=mrate,
                            dom_ok=bool(ok[g])))
        del D
    with open(C.ATT + '/oraclehq_chr%s.json' % ch, 'w') as f:
        json.dump(dict(chrom=str(ch), n_domains=len(doms),
                       wall_s=round(time.time() - t0, 1), per_domain=rep), f)
    with open(C.ATT + '/oraclehq_chr%s.done' % ch, 'w') as f:
        f.write('%d domains %.1f s\n' % (len(doms), time.time() - t0))

if __name__ == '__main__':
    main()
