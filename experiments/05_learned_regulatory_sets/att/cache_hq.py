import json
import os
import sys
import time
import numpy as np
from cyvcf2 import VCF

import att_common as C

R2_MIN = float(os.environ.get('R2_MIN', '0.9'))
OUT = C.FSET + '/att/hqcache'

def main():
    ch = sys.argv[1]
    os.makedirs(OUT, exist_ok=True)
    d = C.master()
    key = np.asarray(d['key'])
    pos19 = np.asarray(d['pos19'])
    r2 = np.asarray(d['r2'], dtype=np.float64)
    chrom = np.asarray(d['chrom'])

    sel = np.where((chrom == str(ch)) & np.isfinite(r2) & (r2 >= R2_MIN))[0]
    sel = sel[np.argsort(pos19[sel])]
    M = len(sel)
    print('chr%s 고품질 변이 %d개' % (ch, M), flush=True)
    if M == 0:
        np.save(OUT + '/hqidx_chr%s.npy' % ch, np.zeros(0, dtype=np.int32))
        open(OUT + '/hq_chr%s.done' % ch, 'w').write('0\n')
        return

    want = {}
    for j, r in enumerate(sel):
        want.setdefault(key[r], j)

    from numpy.lib.format import open_memmap
    out = open_memmap(OUT + '/hq_chr%s.npy' % ch, mode='w+',
                      dtype=np.float32, shape=(M, C.NTOT))
    got = np.zeros(M, dtype=bool)

    sp = np.sort(np.unique(pos19[sel]))
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

    vcf = VCF(C.VCFD + '/chr%s.vcf.gz' % ch)
    assert len(vcf.samples) == C.NTOT, 'sample 수 불일치'
    t0 = time.time()
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
                    ds = np.asarray(v.format('DS'), dtype=np.float32).ravel()
                except Exception:
                    ds = None
                if ds is None or ds.shape[0] != C.NTOT:
                    continue
                bad = ~np.isfinite(ds) | (ds < 0)
                if bad.all():
                    continue
                if bad.any():
                    ds = np.where(bad, np.nanmean(ds[~bad]), ds)
                out[j] = ds
                got[j] = True
    out.flush()
    del out
    np.save(OUT + '/hqidx_chr%s.npy' % ch, sel.astype(np.int32))
    np.save(OUT + '/hqgot_chr%s.npy' % ch, got)
    meta = dict(chrom=str(ch), n_hq=int(M), n_matched=int(got.sum()),
                match_rate=float(got.mean()), r2_min=R2_MIN,
                n_records=int(recs), n_intervals=len(iv),
                wall_s=round(time.time() - t0, 1),
                gb=round(M * C.NTOT * 4 / 1e9, 2))
    json.dump(meta, open(OUT + '/hq_chr%s.json' % ch, 'w'), indent=1)
    open(OUT + '/hq_chr%s.done' % ch, 'w').write(json.dumps(meta) + '\n')
    print('chr%s 완료: 매칭 %d/%d (%.1f%%), %.0f초, %.2f GB'
          % (ch, got.sum(), M, 100 * got.mean(), meta['wall_s'], meta['gb']),
          flush=True)

if __name__ == '__main__':
    main()
