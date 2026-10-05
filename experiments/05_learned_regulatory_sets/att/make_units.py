import collections
import json
import numpy as np

import att_common as C

OUT = C.ATT + '/units.npz'
MIN_VAR = 10
MAX_VAR = 400
SEED = C.SEED + 777

def main():
    d = C.master()
    ti = np.asarray(d['tgt_indices'])
    tp = np.asarray(d['tgt_indptr'])
    ch = np.asarray(d['chrom'])
    p19 = np.asarray(d['pos19'])
    maf = np.asarray(d['maf'], dtype=float)
    nv = len(p19)
    rng = np.random.default_rng(SEED)

    t2v = collections.defaultdict(list)
    for v in range(nv):
        for j in range(tp[v], tp[v + 1]):
            t2v[int(ti[j])].append(v)

    bychr = collections.defaultdict(list)
    for v in range(nv):
        bychr[str(ch[v])].append(v)
    for k in bychr:
        a = np.array(bychr[k])
        bychr[k] = a[np.argsort(p19[a])]

    units = []
    skipped = collections.Counter()
    for t, vs in sorted(t2v.items()):
        vs = np.array(sorted(set(vs)))
        if len(vs) < MIN_VAR:
            skipped['too_few'] += 1
            continue
        cs = set(ch[vs].tolist())
        if len(cs) != 1:
            skipped['multi_chrom'] += 1
            continue
        c = cs.pop()
        if len(vs) > MAX_VAR:
            vs = np.sort(rng.choice(vs, MAX_VAR, replace=False))
            skipped['capped'] += 1
        lo, hi = p19[vs].min(), p19[vs].max()
        pool = bychr[c]
        inwin = pool[(p19[pool] >= lo) & (p19[pool] <= hi)]
        if len(inwin) < len(vs):
            skipped['bad_window'] += 1
            continue
        units.append(('TARGET', int(t), vs))
        w = inwin if len(inwin) <= MAX_VAR else np.sort(
            rng.choice(inwin, MAX_VAR, replace=False))
        units.append(('WINDOW', int(t), w))
        units.append(('NEIGHRND', int(t),
                      np.sort(rng.choice(inwin, len(vs), replace=False))))

    kinds = np.array([u[0] for u in units])
    keys = np.array([u[1] for u in units], dtype=np.int32)
    chs = np.array([str(ch[u[2][0]]) for u in units])
    lens = np.array([len(u[2]) for u in units], dtype=np.int32)
    flat = np.concatenate([u[2] for u in units]).astype(np.int32)
    ptr = np.zeros(len(units) + 1, dtype=np.int64)
    ptr[1:] = np.cumsum(lens)

    print('단위 %d개 (표적 %d개 × 3종)' % (len(units), len(units) // 3))
    print('건너뜀:', dict(skipped))
    for k in ('TARGET', 'WINDOW', 'NEIGHRND'):
        m = kinds == k
        L = lens[m]
        sp = []
        for i in np.where(m)[0][:300]:
            vs = flat[ptr[i]:ptr[i + 1]]
            sp.append((p19[vs].max() - p19[vs].min()) / 1000)
        print('  %-9s n=%4d  변이 중앙 %3d (p90 %3d)  구간 중앙 %.0f kb  '
              'MAF 중앙 %.5f'
              % (k, m.sum(), int(np.median(L)), int(np.percentile(L, 90)),
                 np.median(sp),
                 float(np.median([np.median(maf[flat[ptr[i]:ptr[i + 1]]])
                                  for i in np.where(m)[0][:300]]))))
    cc = collections.Counter()
    for i in range(len(units)):
        cc[chs[i]] += 0
    uniq = collections.defaultdict(set)
    for i in range(len(units)):
        uniq[chs[i]].update(flat[ptr[i]:ptr[i + 1]].tolist())
    mx = max((len(v), k) for k, v in uniq.items())
    print('\n염색체별 고유 변이: 최대 chr%s %d개 -> dosage %.1f GB(float32)'
          % (mx[1], mx[0], mx[0] * C.NTOT * 4 / 1e9))
    print('전 염색체 고유 변이 합 %d' % sum(len(v) for v in uniq.values()))

    np.savez_compressed(OUT, kinds=kinds, keys=keys, chroms=chs,
                        ptr=ptr, flat=flat)
    json.dump(dict(n_units=len(units), n_targets=len(units) // 3,
                   min_var=MIN_VAR, max_var=MAX_VAR, seed=SEED,
                   skipped=dict(skipped)),
              open(C.ATT + '/units_meta.json', 'w'), indent=1,
              ensure_ascii=False)
    print('->', OUT)

if __name__ == '__main__':
    main()
