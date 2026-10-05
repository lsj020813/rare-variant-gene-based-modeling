import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import os
import gzip
import json
import zlib
import numpy as np

FSET = _config_path('${PROJECT_ROOT}/work/fset')
PRIM = FSET + '/primary'
F1F2 = FSET + '/f1f2/ub_inub'
OUT = FSET + '/f3f6'
CACHE = OUT + '/cache'
SEED = 20260922
BOOT = 1000

STAAR_NUM = ['st_cons', 'st_epi_active', 'st_epi_repr', 'st_epi_trans',
             'st_tf', 'st_linsight']
STAAR_TXT = ['st_cage_prom', 'st_genehancer']
GNOM = ['gn_lof_oe', 'gn_loeuf', 'gn_pli', 'gn_lof_z', 'gn_mis_z', 'gn_syn_z']
CADD = ['cadd_raw', 'cadd_phred']
NBIO = 10

def _f(x):
    if x is None:
        return np.nan
    x = x.strip()
    if x == '' or x == 'NA' or x == 'nan' or x == '.' or x == 'None':
        return np.nan
    try:
        return float(x)
    except ValueError:
        return np.nan

def build_cache():
    os.makedirs(CACHE, exist_ok=True)
    key = []
    dom = []
    chrom = []
    pos38 = []
    pos19 = []
    maf = []
    r2 = []
    typed = []
    forced = []
    hccre = []
    ntf = []
    nre2g = []
    dnub = []
    with open(PRIM + '/primary_sample.tsv') as fh:
        h = fh.readline().rstrip('\n').split('\t')
        I = dict((c, i) for i, c in enumerate(h))
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if p[I['in_ub']] != '1':
                continue
            k = p[I['key']]
            key.append(k)
            dom.append(p[I['domain']])
            chrom.append(p[I['chr']])
            pos38.append(int(p[I['pos38']]))
            pos19.append(int(k.split(':')[1]))
            maf.append(_f(p[I['maf']]))
            r2.append(_f(p[I['r2']]))
            typed.append(int(p[I['typed']]))
            forced.append(int(p[I['forced']]))
            hccre.append(_f(p[I['has_ccre']]))
            ntf.append(_f(p[I['n_tf']]))
            nre2g.append(_f(p[I['n_re2g']]))
            dnub.append(int(p[I['domain_n_ub']]))
    n = len(key)
    assert n == 364430, 'in_ub 행 수 불일치: %d' % n
    key = np.array(key)
    dom = np.array(dom)
    d = dict(
        key=key, domain=dom, chrom=np.array(chrom),
        pos38=np.array(pos38, dtype=np.int64),
        pos19=np.array(pos19, dtype=np.int64),
        maf=np.array(maf, dtype=np.float64),
        r2=np.array(r2, dtype=np.float64),
        typed=np.array(typed, dtype=np.int8),
        forced=np.array(forced, dtype=np.int8),
        has_ccre=np.array(hccre, dtype=np.float32),
        n_tf=np.array(ntf, dtype=np.float32),
        n_re2g=np.array(nre2g, dtype=np.float32),
        domain_n_ub=np.array(dnub, dtype=np.int32),
    )

    rid = {}
    for i in range(n):
        rid[(key[i], dom[i])] = i

    combos = ['spectral|eucl', 'ward|eucl', 'average|eucl', 'average|bio',
              'spectral|bio', 'dbscan|bio']
    lab = dict((c, np.full(n, -1, dtype=np.int32)) for c in combos)
    nmatch = dict((c, 0) for c in combos)
    with gzip.open(F1F2 + '/consensus_clusters.csv.gz', 'rt') as fh:
        h = fh.readline().rstrip('\n').split(',')
        J = dict((c, i) for i, c in enumerate(h))
        for line in fh:
            p = line.rstrip('\n').split(',')
            c = p[J['algo']] + '|' + p[J['similarity']]
            if c not in lab:
                continue
            i = rid.get((p[J['key']], p[J['domain']]))
            if i is None:
                continue
            lab[c][i] = int(float(p[J['consensus_cluster']]))
            nmatch[c] += 1
    for c in combos:
        d['mod_' + c.replace('|', '_')] = lab[c]
    assert nmatch['spectral|eucl'] == n, \
        'spectral|eucl 배정 결합 실패: %d/%d' % (nmatch['spectral|eucl'], n)

    vcols = GNOM + STAAR_NUM + CADD
    V = np.full((n, len(vcols)), np.nan, dtype=np.float32)
    TXT = np.zeros((n, len(STAAR_TXT)), dtype=np.int8)
    TXT_seen = np.zeros(n, dtype=np.int8)
    nv = 0
    for ch in [str(x) for x in range(1, 23)]:
        p = PRIM + '/vset_chr%s.tsv.gz' % ch
        if not os.path.exists(p):
            continue
        with gzip.open(p, 'rt') as fh:
            h = fh.readline().rstrip('\n').split('\t')
            K = dict((c, i) for i, c in enumerate(h))
            ci = [K[c] for c in vcols]
            ti = [K[c] for c in STAAR_TXT]
            for line in fh:
                q = line.rstrip('\n').split('\t')
                i = rid.get((q[K['key']], q[K['domain']]))
                if i is None:
                    continue
                nv += 1
                for jj, cc in enumerate(ci):
                    V[i, jj] = _f(q[cc]) if cc < len(q) else np.nan
                TXT_seen[i] = 1
                for jj, cc in enumerate(ti):
                    s = q[cc].strip() if cc < len(q) else ''
                    TXT[i, jj] = 0 if (s == '' or s == 'NA' or s == '.') else 1
    d['vset_cols'] = np.array(vcols)
    d['vset'] = V
    d['staar_txt_cols'] = np.array(STAAR_TXT)
    d['staar_txt'] = TXT
    sidx = [vcols.index(c) for c in STAAR_NUM]
    staar_any = np.isfinite(V[:, sidx]).any(axis=1)
    d['staar_txt_seen'] = (staar_any | (TXT.sum(axis=1) > 0)).astype(np.int8)
    d['vset_row_seen'] = TXT_seen
    d['n_vset_joined'] = np.array([nv])

    SC = np.full((n, NBIO), np.nan, dtype=np.float32)
    DT = np.full((n, NBIO), np.nan, dtype=np.float32)
    vocab = {}
    tgt_rows = [[] for _ in range(n)]
    kidx = {}
    for i in range(n):
        kidx.setdefault(key[i], []).append(i)
    nre = 0
    for ch in [str(x) for x in range(1, 23)]:
        p = PRIM + '/vset_re2g_chr%s.tsv.gz' % ch
        if not os.path.exists(p):
            continue
        with gzip.open(p, 'rt') as fh:
            h = fh.readline().rstrip('\n').split('\t')
            K = dict((c, i) for i, c in enumerate(h))
            sc_i = [K['re2g%d_score' % b] for b in range(NBIO)]
            dt_i = [K['re2g%d_dist_tss' % b] for b in range(NBIO)]
            tg_i = [K['re2g%d_target_ensg' % b] for b in range(NBIO)]
            for line in fh:
                q = line.rstrip('\n').split('\t')
                rr = kidx.get(q[K['key']])
                if rr is None:
                    continue
                nre += 1
                sc = [_f(q[x]) if x < len(q) else np.nan for x in sc_i]
                dt = [_f(q[x]) if x < len(q) else np.nan for x in dt_i]
                tg = set()
                for x in tg_i:
                    if x >= len(q):
                        continue
                    s = q[x].strip()
                    if s == '' or s == 'NA' or s == '.':
                        continue
                    for t in s.split(','):
                        t = t.split('.')[0].strip()
                        if t:
                            if t not in vocab:
                                vocab[t] = len(vocab)
                            tg.add(vocab[t])
                tg = sorted(tg)
                for i in rr:
                    SC[i, :] = sc
                    DT[i, :] = dt
                    tgt_rows[i] = tg
    d['re2g_score'] = SC
    d['re2g_dist_tss'] = DT
    d['n_re2g_joined'] = np.array([nre])
    ind = []
    ptr = [0]
    for i in range(n):
        ind.extend(tgt_rows[i])
        ptr.append(len(ind))
    d['tgt_indices'] = np.array(ind, dtype=np.int32)
    d['tgt_indptr'] = np.array(ptr, dtype=np.int64)
    d['tgt_vocab_n'] = np.array([len(vocab)])

    z = np.load(PRIM + '/Cset_ub_inub.npz', allow_pickle=True)
    zk = z['key']
    zd = z['domain']
    order = np.full(n, -1, dtype=np.int64)
    for j in range(len(zk)):
        i = rid.get((zk[j], zd[j]))
        if i is not None:
            order[i] = j
    assert (order >= 0).all(), 'Cset 행 결합 실패 %d' % int((order < 0).sum())
    d['X'] = np.asarray(z['X'], dtype=np.float32)[order]
    d['X_cols'] = np.asarray(z['columns'])

    udom, dinv = np.unique(dom, return_inverse=True)
    d['udom'] = udom
    d['dom_idx'] = dinv.astype(np.int32)
    ctr = np.zeros(len(udom))
    for g in range(len(udom)):
        m = dinv == g
        ctr[g] = np.median(d['pos38'][m])
    d['dom_center'] = ctr
    dist = np.abs(d['pos38'] - ctr[dinv])
    d['dist_center'] = np.log10(1.0 + dist)
    dch = np.empty(len(udom), dtype=d['chrom'].dtype)
    dfo = np.zeros(len(udom), dtype=np.int8)
    dnu = np.zeros(len(udom), dtype=np.int32)
    for g in range(len(udom)):
        m = np.where(dinv == g)[0]
        dch[g] = d['chrom'][m[0]]
        dfo[g] = d['forced'][m[0]]
        dnu[g] = len(m)
    d['dom_chrom'] = dch
    d['dom_forced'] = dfo
    d['dom_n'] = dnu

    np.savez_compressed(CACHE + '/master.npz', **d)
    meta = dict(n_rows=int(n), n_domains=int(len(udom)),
                n_vset_joined=int(nv), n_re2g_joined=int(nre),
                tgt_vocab=int(len(vocab)),
                mod_match=dict((k, int(v)) for k, v in nmatch.items()))
    with open(CACHE + '/master_meta.json', 'w') as fh:
        json.dump(meta, fh, indent=1)
    return meta

_M = {}

def load():
    if 'd' not in _M:
        z = np.load(CACHE + '/master.npz', allow_pickle=False)
        _M['d'] = dict((k, z[k]) for k in z.files)
    return _M['d']

CHR_LEN = {'1': 248956422, '2': 242193529, '3': 198295559, '4': 190214555,
           '5': 181538259, '6': 170805979, '7': 159345973, '8': 145138636,
           '9': 138394717, '10': 133797422, '11': 135086622, '12': 133275309,
           '13': 114364328, '14': 107043718, '15': 101991189, '16': 90338345,
           '17': 83257441, '18': 80373285, '19': 58617616, '20': 64444167,
           '21': 46709983, '22': 50818468}

def domain_strata(d):
    ch = [c.replace('chr', '') for c in d['dom_chrom']]
    L = np.array([CHR_LEN.get(c, 1) for c in ch], dtype=float)
    ct = np.quantile(L, [1 / 3.0, 2 / 3.0])
    a = np.digitize(L, ct)
    nq = np.quantile(d['dom_n'], [0.25, 0.5, 0.75])
    b = np.digitize(d['dom_n'], nq)
    return a * 4 + b

def pick_domains(d, n_pick, random_only=True, seed=SEED):
    st = domain_strata(d)
    elig = np.where(d['dom_forced'] == 0)[0] if random_only \
        else np.arange(len(d['udom']))
    rng = np.random.default_rng(seed)
    sel = []
    us, cnt = np.unique(st[elig], return_counts=True)
    quota = np.maximum(1, np.round(cnt / float(cnt.sum()) * n_pick).astype(int))
    for s, q in zip(us, quota):
        pool = elig[st[elig] == s]
        q = int(min(q, len(pool)))
        sel.extend(rng.choice(pool, q, replace=False).tolist())
    sel = np.array(sorted(set(sel)))
    if len(sel) > n_pick:
        sel = np.array(sorted(rng.choice(sel, n_pick, replace=False).tolist()))
    return sel

def _terc(v):
    ok = np.isfinite(v)
    if ok.sum() < 3:
        return np.zeros(len(v), dtype=np.int64)
    q = np.quantile(v[ok], [1 / 3.0, 2 / 3.0])
    return np.digitize(np.where(ok, v, np.nanmedian(v[ok])), q)

def _med(v):
    ok = np.isfinite(v)
    if ok.sum() == 0:
        return np.zeros(len(v), dtype=np.int64)
    m = np.median(v[ok])
    return (np.where(ok, v, m) > m).astype(np.int64)

def strata_for(d, rows, kind, ld_ctx=None):
    re2 = (d['n_re2g'][rows] > 0).astype(np.int64)
    if kind == 'A':
        return np.zeros(len(rows), dtype=np.int64)
    if kind == 'B':
        return (((re2 * 3 + _terc(d['maf'][rows])) * 3
                 + _terc(d['r2'][rows])) * 3 + _terc(d['dist_center'][rows]))
    if kind == 'C':
        cc = (d['has_ccre'][rows] > 0).astype(np.int64)
        tf = (d['n_tf'][rows] > 0).astype(np.int64)
        return (((((re2 * 2 + _med(d['maf'][rows])) * 2
                   + _med(d['r2'][rows])) * 2
                  + _med(d['dist_center'][rows])) * 2 + cc) * 2 + tf)
    if kind == 'M':
        v = d['maf'][rows]
        ok = np.isfinite(v)
        if ok.sum() >= 10:
            q = np.quantile(v[ok], np.arange(1, 10) / 10.0)
            b = np.digitize(np.where(ok, v, np.median(v[ok])), q)
        else:
            b = np.zeros(len(rows), dtype=np.int64)
        return re2 * 10 + b
    if kind == 'LD':
        assert ld_ctx is not None
        return (((re2 * 3 + _terc(d['maf'][rows])) * 3
                 + _terc(d['dist_center'][rows])) * 3 + _terc(ld_ctx))
    raise ValueError(kind)

def det_seed(*parts):
    s = '|'.join([str(p) for p in parts]).encode('utf-8')
    return (zlib.crc32(s) ^ SEED) % (2 ** 31 - 1)

def matched_null(pool_strata, member_mask, boot=BOOT, seed=0):
    rng = np.random.default_rng(seed)
    nm = int(member_mask.sum())
    out = np.empty((boot, nm), dtype=np.int64)
    off = 0
    for c in np.unique(pool_strata[member_mask]):
        k = int((pool_strata[member_mask] == c).sum())
        cd = np.where(pool_strata == c)[0]
        p = len(cd)
        if p <= k:
            out[:, off:off + k] = cd[None, :][:, :k].repeat(boot, axis=0)
        elif k == 1:
            out[:, off:off + 1] = cd[rng.integers(0, p, size=(boot, 1))]
        else:
            r = rng.random((boot, p))
            sel = np.argpartition(r, k - 1, axis=1)[:, :k]
            out[:, off:off + k] = cd[sel]
        off += k
    assert off == nm
    return out
