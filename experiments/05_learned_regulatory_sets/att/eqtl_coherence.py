import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import collections
import csv
import json
import sys
import numpy as np

sys.path.insert(0, _config_path('${PROJECT_ROOT}/work/fset/att'))
import att_common as C

NRAND = 20
MINV = 8

def main():
    eq = sys.argv[1]
    tag = sys.argv[2]
    d = C.master()
    key = np.asarray(d['key'])
    p38 = np.asarray(d['pos38'])
    ch = np.asarray(d['chrom'])
    di = np.asarray(d['dom_idx'])
    ud = np.asarray(d['udom'])
    hcc = np.asarray(d['has_ccre'], dtype=float) > 0
    maf = np.asarray(d['maf'], dtype=float)

    gene_of_dom = [str(x).split('.')[0] for x in ud]
    alle = {}
    for v in range(len(key)):
        parts = str(key[v]).split(':')
        alle[v] = (parts[2], parts[3])
    loc = collections.defaultdict(list)
    for v in range(len(key)):
        loc[(str(ch[v]), int(p38[v]))].append(v)

    zmap = collections.defaultdict(dict)
    n_rows = n_hit = 0
    with open(eq) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 17:
                continue
            n_rows += 1
            g = f[16]
            cand = loc.get((f[1], int(f[2])))
            if not cand:
                continue
            ref, alt = f[3], f[4]
            try:
                b = float(f[9]); se = float(f[10])
            except ValueError:
                continue
            if not np.isfinite(b) or not np.isfinite(se) or se <= 0:
                continue
            for v in cand:
                r0, a0 = alle[v]
                if (r0, a0) == (ref, alt):
                    zmap[g][v] = b / se
                elif (r0, a0) == (alt, ref):
                    zmap[g][v] = -b / se
                else:
                    continue
                n_hit += 1
    print('eQTL 행 %d 처리, 매칭 %d, 유전자 %d개' % (n_rows, n_hit, len(zmap)),
          flush=True)

    rng = np.random.default_rng(C.SEED + 2718)
    rows = []
    for g in range(C.NDOM):
        gene = gene_of_dom[g]
        zz = zmap.get(gene)
        if not zz:
            continue
        rid = np.where(di == g)[0]
        have = np.array([v for v in rid if v in zz])
        if len(have) < 2 * MINV:
            continue
        z = np.array([zz[v] for v in have])
        cc = hcc[have]
        for nm, m in (('CCRE', cc), ('NONCCRE', ~cc)):
            if m.sum() < MINV:
                continue
            zs = z[m]
            rows.append(dict(domain=g, gene=gene, unit=nm, n=int(m.sum()),
                             mean_abs_z=float(np.mean(np.abs(zs))),
                             coherence=float(abs(zs.sum())
                                              / np.sqrt(len(zs)
                                                        * (zs * zs).sum())),
                             maf_median=float(np.median(maf[have][m]))))
            k = int(m.sum())
            for _ in range(NRAND):
                idx = rng.choice(len(z), k, replace=False)
                zs2 = z[idx]
                rows.append(dict(domain=g, gene=gene, unit='RAND_' + nm,
                                 n=k,
                                 mean_abs_z=float(np.mean(np.abs(zs2))),
                                 coherence=float(abs(zs2.sum())
                                                 / np.sqrt(len(zs2)
                                                           * (zs2 * zs2).sum())),
                                 maf_median=float(np.median(maf[have][idx]))))
    print('단위 레코드 %d개 (도메인 %d개)'
          % (len(rows), len(set(r['domain'] for r in rows))), flush=True)
    if not rows:
        return

    out = C.ATT + '/eqtl_%s' % tag
    with open(out + '_units.csv', 'w') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)

    res = []
    for nm in ('CCRE', 'RAND_CCRE', 'NONCCRE', 'RAND_NONCCRE'):
        s = [r for r in rows if r['unit'] == nm]
        if not s:
            continue
        a = np.array([r['mean_abs_z'] for r in s])
        cvals = np.array([r['coherence'] for r in s])
        res.append(dict(unit=nm, n_rec=len(s),
                        n_var_median=float(np.median([r['n'] for r in s])),
                        mean_abs_z=float(np.median(a)),
                        coherence=float(np.median(cvals))))
    h = '%-14s %8s %10s %13s %12s'
    print()
    print(h % ('단위', 'n레코드', '변이중앙', '평균|z| 중앙', '정합성 중앙'))
    for r in res:
        print(h % (r['unit'], r['n_rec'], '%.0f' % r['n_var_median'],
                   '%.4f' % r['mean_abs_z'], '%.4f' % r['coherence']))
    for base in ('CCRE', 'NONCCRE'):
        byd = collections.defaultdict(dict)
        for r in rows:
            if r['unit'] == base:
                byd[r['domain']]['obs'] = r
            elif r['unit'] == 'RAND_' + base:
                byd[r['domain']].setdefault('rnd', []).append(r)
        da, dc = [], []
        for g, v in byd.items():
            if 'obs' not in v or 'rnd' not in v:
                continue
            da.append(v['obs']['mean_abs_z']
                      - np.mean([x['mean_abs_z'] for x in v['rnd']]))
            dc.append(v['obs']['coherence']
                      - np.mean([x['coherence'] for x in v['rnd']]))
        if da:
            da = np.array(da); dc = np.array(dc)
            print()
            print('%s − 같은도메인 무작위 (도메인 %d개):' % (base, len(da)))
            print('  Δ평균|z|  중앙 %+.4f  >0 비율 %.1f%%'
                  % (np.median(da), 100 * np.mean(da > 0)))
            print('  Δ정합성   중앙 %+.4f  >0 비율 %.1f%%'
                  % (np.median(dc), 100 * np.mean(dc > 0)))
            res.append(dict(unit=base + '_vs_random',
                            d_mean_abs_z=float(np.median(da)),
                            d_coherence=float(np.median(dc)),
                            frac_dz_pos=float(np.mean(da > 0)),
                            frac_dc_pos=float(np.mean(dc > 0)),
                            n_dom=len(da)))
    json.dump(dict(tag=tag, n_rows=n_rows, n_hit=n_hit, result=res),
              open(out + '_summary.json', 'w'), indent=1, ensure_ascii=False)
    print('\n->', out + '_summary.json')

if __name__ == '__main__':
    main()
