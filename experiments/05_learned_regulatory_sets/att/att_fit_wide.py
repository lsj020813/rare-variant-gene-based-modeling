import sys
import csv
import json
import time
import numpy as np
import torch
from scipy.stats import norm

import att_common as C

import os as _os
OUT = C.ATT + '/wide'
PRIVOUT = OUT + '/private'
for _p in (OUT, PRIVOUT):
    _os.makedirs(_p, exist_ok=True)

torch.set_num_threads(4)
NFOLD = 5
NINNER = 5
LAM = [0.0, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5]
EPOCHS = 400
LR = 0.02
HID = 16
WLO, WHI = -0.5, 1.5
NPC = 10

def person_folds():
    rng = np.random.default_rng(C.SEED)
    p = rng.permutation(C.NTOT)
    f = np.empty(C.NTOT, dtype=np.int8)
    for k in range(NFOLD):
        f[p[k::NFOLD]] = k
    return f

def load_pheno():
    from cyvcf2 import VCF
    v = VCF(C.VCFD + '/chr22.vcf.gz')
    pos = dict((s, i) for i, s in enumerate(v.samples))
    idx, ys, cs = [], [], []
    with open(C.PHENO) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            i = pos.get(r['sample_id'])
            if i is None:
                continue
            try:
                yv = float(r['y'])
                cv = [float(r[c]) for c in C.COVARS]
            except (ValueError, TypeError, KeyError):
                continue
            if not np.isfinite(yv) or not np.isfinite(cv).all():
                continue
            idx.append(i)
            ys.append(yv)
            cs.append(cv)
    idx = np.asarray(idx, dtype=np.int64)
    o = np.argsort(idx)
    return (idx[o], np.asarray(ys, dtype=np.float64)[o],
            np.asarray(cs, dtype=np.float64)[o])

def feature_matrix(raw, cmean):
    R = raw[:, :, :5].reshape(C.NDOM * 2, 5).astype(np.float64)
    Cm = cmean.reshape(C.NDOM * 2, -1).astype(np.float64)
    m = np.isfinite(R).all(1) & np.isfinite(Cm).all(1)
    Xc = Cm[m]
    sd = Xc.std(0)
    sd[sd < 1e-12] = 1.0
    Zc = (Xc - Xc.mean(0)) / sd
    _, _, Vt = np.linalg.svd(Zc, full_matrices=False)
    F = np.zeros((C.NDOM * 2, 5 + NPC))
    F[m, :5] = R[m]
    F[m, 5:] = Zc @ Vt[:NPC].T
    sd2 = F[m].std(0)
    sd2[sd2 < 1e-12] = 1.0
    F = (F - F[m].mean(0)) / sd2
    F[~m] = 0.0
    return F.reshape(C.NDOM, 2, 5 + NPC), m.reshape(C.NDOM, 2).all(1)

def make_net(d, hid, seed):
    torch.manual_seed(seed)
    if hid:
        return torch.nn.Sequential(torch.nn.Linear(d, hid), torch.nn.ReLU(),
                                   torch.nn.Linear(hid, 1))
    return torch.nn.Linear(d, 1)

def gate_w(net, Xmaj, Xmin):
    s = torch.sigmoid(net(Xmaj).squeeze(-1) - net(Xmin).squeeze(-1))
    return WLO + (WHI - WLO) * s

def ev_terms(w, st):
    u, v, r, va, vb = st
    num = w * u + (1.0 - w) * v
    den = w * w * va + (1.0 - w) ** 2 * vb + 2.0 * w * (1.0 - w) * r
    return num * num / torch.clamp(den, min=1e-6)

def fit_gate(Xmj, Xmn, st, msk, lam, hid, seed):
    net = make_net(Xmj.shape[1], hid, seed)
    opt = torch.optim.Adam(net.parameters(), lr=LR)
    ng = float(msk.sum())
    for _ in range(EPOCHS):
        opt.zero_grad()
        ev = ev_terms(gate_w(net, Xmj, Xmn), st)
        reg = sum((p * p).sum() for p in net.parameters())
        (-(ev * msk).sum() / ng + lam * reg).backward()
        opt.step()
    return net

def eval_gate(net, Xmj, Xmn, st, msk):
    with torch.no_grad():
        return float((ev_terms(gate_w(net, Xmj, Xmn), st) * msk).sum()
                     / float(msk.sum()))

def block_stats(A, B, yt, blocks):
    nb = len(blocks)
    G = A.shape[0]
    o = dict(n=np.zeros(nb), Sy=np.zeros(nb),
             Sa=np.zeros((nb, G)), Saa=np.zeros((nb, G)),
             Say=np.zeros((nb, G)), Sb=np.zeros((nb, G)),
             Sbb=np.zeros((nb, G)), Sby=np.zeros((nb, G)),
             Sab=np.zeros((nb, G)))
    for k, cols in enumerate(blocks):
        a = A[:, cols].astype(np.float64)
        b = B[:, cols].astype(np.float64)
        yk = yt[cols]
        o['n'][k] = len(cols)
        o['Sy'][k] = yk.sum()
        o['Sa'][k] = a.sum(1)
        o['Saa'][k] = np.einsum('ij,ij->i', a, a)
        o['Say'][k] = a @ yk
        o['Sb'][k] = b.sum(1)
        o['Sbb'][k] = np.einsum('ij,ij->i', b, b)
        o['Sby'][k] = b @ yk
        o['Sab'][k] = np.einsum('ij,ij->i', a, b)
        del a, b
    return o

def moments(o, ks):
    n = o['n'][ks].sum()
    my = o['Sy'][ks].sum() / n
    ma = o['Sa'][ks].sum(0) / n
    mb = o['Sb'][ks].sum(0) / n
    va = o['Saa'][ks].sum(0) / n - ma * ma
    vb = o['Sbb'][ks].sum(0) / n - mb * mb
    cay = o['Say'][ks].sum(0) / n - ma * my
    cby = o['Sby'][ks].sum(0) / n - mb * my
    cab = o['Sab'][ks].sum(0) / n - ma * mb
    return dict(n=n, ma=ma, mb=mb, va=np.maximum(va, 0.0),
                vb=np.maximum(vb, 0.0), cay=cay, cby=cby, cab=cab)

def to_st(m, sa, sb, oa, ob):
    da = np.where(oa, sa, 1.0)
    db = np.where(ob, sb, 1.0)
    u = np.where(oa, m['cay'] / da, 0.0)
    v = np.where(ob, m['cby'] / db, 0.0)
    r = np.where(oa & ob, m['cab'] / (da * db), 0.0)
    va = np.where(oa, m['va'] / (da * da), 0.0)
    vb = np.where(ob, m['vb'] / (db * db), 0.0)
    return tuple(torch.tensor(x, dtype=torch.float32)
                 for x in (u, v, r, va, vb))

def resid_y(y, Xc, tr):
    D = np.column_stack([np.ones(len(y)), Xc])
    beta, *_ = np.linalg.lstsq(D[tr], y[tr], rcond=None)
    e = y - D @ beta
    s = e[tr].std()
    e = (e - e[tr].mean()) / (s if s > 1e-12 else 1.0)
    return e

def std_apply(A, cols, m, s, ok):
    a = A[:, cols]
    return np.where(ok[:, None],
                    (a - m[:, None]) / np.where(ok, s, 1.0)[:, None], 0.0)

def zstats(S, y, Xc):
    D = np.column_stack([np.ones(len(y)), Xc])
    Q, _ = np.linalg.qr(D)
    yp = y - Q @ (Q.T @ y)
    Sp = S - (S @ Q) @ Q.T
    ss = np.einsum('ij,ij->i', Sp, Sp)
    good = ss > 1e-8
    num = Sp @ yp
    beta = np.where(good, num / np.where(good, ss, 1.0), np.nan)
    dfree = len(y) - D.shape[1] - 1
    rss = float(yp @ yp) - np.where(good, num * num / np.where(good, ss, 1.0),
                                    0.0)
    se = np.sqrt(np.where(good, rss / dfree / np.where(good, ss, 1.0), np.nan))
    return beta, beta / se

def bh(p):
    p = np.asarray(p, dtype=np.float64)
    m = np.isfinite(p)
    q = np.full(len(p), np.nan)
    pp = p[m]
    n = len(pp)
    if n == 0:
        return q
    qq = pp * n / (np.argsort(np.argsort(pp)) + 1.0)
    order = np.argsort(-pp)
    run = np.inf
    out = np.empty(n)
    for i in order:
        run = min(run, qq[i])
        out[i] = run
    q[m] = np.minimum(out, 1.0)
    return q

def run_att(Amaj, Amin, F, gmask, y, Xc, fold, lam_fix=None, gfold=None):
    N = len(y)
    Xmj = torch.tensor(F[:, 0, :], dtype=torch.float32)
    Xmn = torch.tensor(F[:, 1, :], dtype=torch.float32)
    S = np.zeros((C.NDOM, N), dtype=np.float32)
    Wacc = np.zeros(C.NDOM)
    sel = []
    for f in range(NFOLD):
        tr = np.where(fold != f)[0]
        te = np.where(fold == f)[0]
        yt = resid_y(y, Xc, tr)
        rng = np.random.default_rng(C.SEED + 1000 + f)
        ip = rng.permutation(len(tr))
        blocks = [tr[ip[k::NINNER]] for k in range(NINNER)]
        o = block_stats(Amaj, Amin, yt, blocks)
        mfull = moments(o, list(range(NINNER)))
        sa = np.sqrt(mfull['va'])
        sb = np.sqrt(mfull['vb'])
        oa = sa > 1e-9
        ob = sb > 1e-9
        stf = to_st(mfull, sa, sb, oa, ob)
        base = gmask & oa & ob
        if lam_fix is None:
            res = {}
            for hid in (0, HID):
                for lam in LAM:
                    sc = []
                    for k in range(NINNER):
                        ks = [j for j in range(NINNER) if j != k]
                        mi = moments(o, ks)
                        sai = np.sqrt(mi['va'])
                        sbi = np.sqrt(mi['vb'])
                        oai = sai > 1e-9
                        obi = sbi > 1e-9
                        mk = torch.tensor((gmask & oai & obi).astype(
                            np.float32))
                        net = fit_gate(Xmj, Xmn,
                                       to_st(mi, sai, sbi, oai, obi), mk,
                                       lam, hid, C.SEED + f * 97 + k)
                        mv = moments(o, [k])
                        sc.append(eval_gate(net, Xmj, Xmn,
                                            to_st(mv, sai, sbi, oai, obi),
                                            mk))
                    res[(hid, lam)] = (float(np.mean(sc)),
                                       float(np.std(sc, ddof=1)
                                             / np.sqrt(NINNER)))
            best = {}
            for hid in (0, HID):
                bm = max(res[(hid, l)][0] for l in LAM)
                bl = [l for l in LAM if res[(hid, l)][0] == bm][0]
                bse = res[(hid, bl)][1]
                best[hid] = (bm, bl, bse)
            lin_m, lin_l, _ = best[0]
            mlp_m, mlp_l, mlp_se = best[HID]
            hid_use, lam_use = (HID, mlp_l)
            sel.append(dict(fold=f, hid=int(hid_use), lam=float(lam_use),
                            inner_linear=lin_m, inner_mlp=mlp_m,
                            inner_linear_lam=lin_l, mlp_se=mlp_se,
                            rule='CV-min on inner 5-fold; MLP 게이트 강제 (A5-CVMIN 변형)'))
        else:
            hid_use, lam_use = lam_fix[f]
            sel.append(dict(fold=f, hid=int(hid_use), lam=float(lam_use),
                            rule='outer λ 재사용'))
        if gfold is None:
            net = fit_gate(Xmj, Xmn, stf,
                           torch.tensor(base.astype(np.float32)),
                           lam_use, hid_use, C.SEED + f)
            with torch.no_grad():
                w = gate_w(net, Xmj, Xmn).numpy().astype(np.float64)
            Wacc += w
            S[:, te] = (w[:, None] * std_apply(Amaj, te, mfull['ma'], sa, oa)
                        + (1 - w)[:, None] * std_apply(Amin, te, mfull['mb'],
                                                       sb, ob))
        else:
            Da = std_apply(Amaj, te, mfull['ma'], sa, oa)
            Db = std_apply(Amin, te, mfull['mb'], sb, ob)
            for h in range(NFOLD):
                net = fit_gate(Xmj, Xmn, stf,
                               torch.tensor((base & (gfold != h)).astype(
                                   np.float32)), lam_use, hid_use,
                               C.SEED + f * 13 + h)
                with torch.no_grad():
                    w = gate_w(net, Xmj, Xmn).numpy().astype(np.float64)
                gi = np.where(gfold == h)[0]
                S[np.ix_(gi, te)] = (w[gi, None] * Da[gi]
                                     + (1 - w[gi, None]) * Db[gi])
        del o
    return S, Wacc / NFOLD, sel

def arm_single(A, y, fold):
    S = np.zeros((C.NDOM, len(y)), dtype=np.float32)
    for f in range(NFOLD):
        tr = np.where(fold != f)[0]
        te = np.where(fold == f)[0]
        a = A[:, tr]
        m = a.mean(1)
        s = a.std(1)
        S[:, te] = std_apply(A, te, m, s, s > 1e-9)
        del a
    return S

def arm_fixed_w(Amaj, Amin, w, y, fold):
    S = np.zeros((C.NDOM, len(y)), dtype=np.float32)
    for f in range(NFOLD):
        tr = np.where(fold != f)[0]
        te = np.where(fold == f)[0]
        a = Amaj[:, tr]
        b = Amin[:, tr]
        ma, sa = a.mean(1), a.std(1)
        mb, sb = b.mean(1), b.std(1)
        del a, b
        S[:, te] = (w[:, None] * std_apply(Amaj, te, ma, sa, sa > 1e-9)
                    + (1 - w)[:, None] * std_apply(Amin, te, mb, sb,
                                                   sb > 1e-9))
    return S

def _int_cols(tot, Amaj, Amin, idx):
    a = np.asarray(Amaj[:, idx], dtype=np.float32)
    b = np.asarray(Amin[:, idx], dtype=np.float32)
    return [np.asarray(tot[:, idx], dtype=np.float32), a * b, a * a, b * b]

def _std_cols(Z, mu, sd, ok):
    return [np.where(o[:, None], (z - m[:, None])
                     / np.where(o, s, 1.0)[:, None], 0.0).astype(np.float32)
            for z, m, s, o in zip(Z, mu, sd, ok)]

def arm_int(tot, Amaj, Amin, y, Xc, fold):
    N = len(y)
    S = np.zeros((C.NDOM, N), dtype=np.float32)
    SP = np.zeros((C.NDOM, N), dtype=np.float32)
    gi = np.arange(C.NDOM)
    ix = [0, 2, 3]
    for f in range(NFOLD):
        tr = np.where(fold != f)[0]
        te = np.where(fold == f)[0]
        yt = resid_y(y, Xc, tr)
        Z = _int_cols(tot, Amaj, Amin, tr)
        mu = [z.mean(1).astype(np.float64) for z in Z]
        sd = [z.std(1).astype(np.float64) for z in Z]
        ok = [s > 1e-9 for s in sd]
        Z = _std_cols(Z, mu, sd, ok)
        ytr = (yt[tr] - yt[tr].mean()).astype(np.float32)
        n = float(len(tr))
        Gm = np.zeros((C.NDOM, 4, 4))
        Xty = np.zeros((C.NDOM, 4))
        for i in range(4):
            Xty[:, i] = (Z[i] @ ytr) / n
            for j in range(i, 4):
                v = np.einsum('gn,gn->g', Z[i], Z[j]) / n
                Gm[:, i, j] = v
                Gm[:, j, i] = v
        beta = np.einsum('gij,gj->gi', np.linalg.pinv(Gm), Xty)
        gam = np.einsum('gij,gj->gi',
                        np.linalg.pinv(Gm[np.ix_(gi, ix, ix)]), Gm[:, ix, 1])
        del Z
        Zt = _std_cols(_int_cols(tot, Amaj, Amin, te), mu, sd, ok)
        acc = np.zeros((C.NDOM, len(te)), dtype=np.float32)
        for i in range(4):
            acc += beta[:, i:i + 1].astype(np.float32) * Zt[i]
        S[:, te] = acc
        pp = Zt[1].copy()
        for k, i in enumerate(ix):
            pp -= gam[:, k:k + 1].astype(np.float32) * Zt[i]
        SP[:, te] = pp
        del Zt, acc, pp
    return S, SP

def load_col(i, cols):
    return np.asarray(np.load(C.BURD + '/col%02d.npy' % i,
                              mmap_mode='r')[:, cols], dtype=np.float32)

def _f(x):
    try:
        x = float(x)
    except (TypeError, ValueError):
        return 'NA'
    return 'NA' if not np.isfinite(x) else repr(round(x, 6))

def main():
    nshuf = int(sys.argv[1]) if len(sys.argv) > 1 else C.NSHUF
    t0 = time.time()
    C.ensure_dirs()
    d = C.master()
    udom = np.asarray(d['udom'])
    dom_chrom = np.asarray(d['dom_chrom'])
    dom_forced = np.asarray(d['dom_forced'])
    dn = np.asarray(d['dom_n'])
    di = np.asarray(d['dom_idx'])
    hcc = np.asarray(d['has_ccre']) > 0
    nccre = np.bincount(di[hcc], minlength=C.NDOM)
    Z = np.load(C.ATT + '/att_labels.npz', allow_pickle=True)
    RAW, CM, NK, gok = Z['raw'], Z['cmean'], Z['nk'], Z['dom_ok']
    fold_all = person_folds()
    np.save(PRIVOUT + '/person_fold.npy', fold_all)
    idx, y, Xc = load_pheno()
    fold = fold_all[idx]
    N = len(y)
    print('[pheno] N=%d fold counts=%s' % (N, np.bincount(fold).tolist()),
          flush=True)
    totlab = load_col(C.C_TOT_LAB, idx)
    maj = load_col(C.C_MAJ0, idx)
    mino = totlab - maj
    F0, gm0 = feature_matrix(RAW[0], CM[0])
    gmask = gm0 & gok
    S_att, w_att, sel = run_att(maj, mino, F0, gmask, y, Xc, fold)
    lam_fix = [(s['hid'], s['lam']) for s in sorted(sel, key=lambda r: r['fold'])]
    print('[obs] ATT done %.1fs sel=%s' % (time.time() - t0, lam_fix),
          flush=True)
    out = {}
    tot = load_col(C.C_TOT_ALL, idx)
    out['FLAT'] = zstats(arm_single(tot, y, fold), y, Xc)
    ccre = load_col(C.C_CCRE, idx)
    out['FIXED'] = zstats(arm_single(ccre, y, fold), y, Xc)
    del ccre
    tt = NK.sum(1).astype(np.float64)
    w_size = np.where(tt > 0, NK[:, 0] / np.where(tt > 0, tt, 1), 0.5)
    out['SIZE'] = zstats(arm_fixed_w(maj, mino, w_size, y, fold), y, Xc)
    out['ATT'] = zstats(S_att, y, Xc)
    _one = np.ones(C.NDOM)
    _zero = np.zeros(C.NDOM)
    S_maj = arm_fixed_w(maj, mino, _one, y, fold)
    S_min = arm_fixed_w(maj, mino, _zero, y, fold)
    out['MAJONLY'] = zstats(S_maj, y, Xc)
    out['MINONLY'] = zstats(S_min, y, Xc)
    out['CONTRAST'] = zstats(S_maj - S_min, y, Xc)
    del S_maj, S_min
    del S_att
    S_int, S_prod = arm_int(tot, maj, mino, y, Xc, fold)
    out['INT'] = zstats(S_int, y, Xc)
    out['INT_product_only'] = zstats(S_prod, y, Xc)
    del S_int, S_prod
    print('[obs] INT done %.1fs' % (time.time() - t0), flush=True)
    gfold = np.argsort(np.random.default_rng(C.SEED).permutation(C.NDOM)
                       ).astype(np.int64) % NFOLD
    S_gh, _, _ = run_att(maj, mino, F0, gmask, y, Xc, fold, lam_fix=lam_fix,
                         gfold=gfold)
    out['ATT_GENEHOLD'] = zstats(S_gh, y, Xc)
    del S_gh
    print('[obs] gene-holdout done %.1fs' % (time.time() - t0), flush=True)
    tl1 = load_col(C.C_TOT001_LAB, idx)
    m1 = load_col(C.C_MAJ001, idx)
    n1 = tl1 - m1
    S1, w_a1, _ = run_att(m1, n1, F0, gmask, y, Xc, fold, lam_fix=lam_fix)
    out['ATT_MAF001'] = zstats(S1, y, Xc)
    del S1, tl1
    t1 = load_col(C.C_TOT001_ALL, idx)
    out['FLAT_MAF001'] = zstats(arm_single(t1, y, fold), y, Xc)
    del t1
    c1 = load_col(C.C_CCRE001, idx)
    out['FIXED_MAF001'] = zstats(arm_single(c1, y, fold), y, Xc)
    del c1, m1, n1
    t_obs = time.time() - t0
    print('[obs] all arms %.1fs' % t_obs, flush=True)
    zflat = out['FLAT'][1]
    dz_obs = np.abs(out['ATT'][1]) - np.abs(zflat)
    rows = []
    for L in range(1, nshuf + 1):
        ts = time.time()
        mj = load_col(C.C_MAJ0 + L, idx)
        mn = totlab - mj
        FL, gmL = feature_matrix(RAW[L], CM[L])
        SL, wL, selL = run_att(mj, mn, FL, gmL & gok, y, Xc, fold)
        bL, zL = zstats(SL, y, Xc)
        dL = np.abs(zL) - np.abs(zflat)
        qL = bh(2.0 * norm.sf(np.abs(zL)))
        _o = np.ones(C.NDOM)
        _z = np.zeros(C.NDOM)
        _sm = arm_fixed_w(mj, mn, _o, y, fold)
        _sn = arm_fixed_w(mj, mn, _z, y, fold)
        _, zCL = zstats(_sm - _sn, y, Xc)
        qCL = bh(2.0 * norm.sf(np.abs(zCL)))
        dCL = np.abs(zCL) - np.abs(zflat)
        del _sm, _sn
        SIL, SPL = arm_int(tot, mj, mn, y, Xc, fold)
        _, zIL = zstats(SIL, y, Xc)
        dIL = np.abs(zIL) - np.abs(zflat)
        qIL = bh(2.0 * norm.sf(np.abs(zIL)))
        del SIL, SPL
        rows.append(dict(shuffle=L, seed=C.SEED + L,
                         mean_delta=float(np.nanmean(dL)),
                         median_delta=float(np.nanmedian(dL)),
                         frac_delta_pos=float(np.nanmean(dL > 0)),
                         n_sig_fdr05=int(np.nansum(qL < 0.05)),
                         max_abs_z=float(np.nanmax(np.abs(zL))),
                         mean_w_major=float(np.nanmean(wL)),
                         hid=int(selL[0]['hid']), lam=float(selL[0]['lam']),
                         mean_delta_int=float(np.nanmean(dIL)),
                         median_delta_int=float(np.nanmedian(dIL)),
                         frac_delta_int_pos=float(np.nanmean(dIL > 0)),
                         n_sig_int_fdr05=int(np.nansum(qIL < 0.05)),
                         n_sig_contrast_fdr05=int(np.nansum(qCL < 0.05)),
                         mean_delta_contrast=float(np.nanmean(dCL)),
                         max_abs_z_contrast=float(np.nanmax(np.abs(zCL))),
                         max_abs_z_int=float(np.nanmax(np.abs(zIL))),
                         wall_s=round(time.time() - ts, 1)))
        np.save(PRIVOUT + '/null_delta_int_%02d.npy' % L,
                dIL.astype(np.float32))
        np.save(PRIVOUT + '/null_delta_%02d.npy' % L, dL.astype(np.float32))
        np.save(PRIVOUT + '/null_z_%02d.npy' % L, zL.astype(np.float32))
        del mj, mn, SL
        print('[null] L=%d %.1fs mean_delta=%.4f nsig=%d'
              % (L, time.time() - ts, rows[-1]['mean_delta'],
                 rows[-1]['n_sig_fdr05']), flush=True)
    np.savez(PRIVOUT + '/att_raw.npz',
             **dict([('z_' + k, v[1]) for k, v in out.items()]),
             **dict([('b_' + k, v[0]) for k, v in out.items()]),
             w_att=w_att, w_size=w_size, w_a1=w_a1, dz_obs=dz_obs,
             dz_int_obs=np.abs(out['INT'][1]) - np.abs(zflat),
             gfold=gfold, gmask=gmask, nk=NK, raw0=RAW[0], nccre=nccre)
    with open(OUT + '/att_label_shuffle_null.csv', 'w') as fh:
        wtr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wtr.writeheader()
        for r in rows:
            wtr.writerow(r)
    arms = ['FLAT', 'FIXED', 'SIZE', 'ATT', 'INT', 'INT_product_only',
            'ATT_GENEHOLD', 'ATT_MAF001', 'FLAT_MAF001', 'FIXED_MAF001',
            'MAJONLY', 'MINONLY', 'CONTRAST']
    qs = dict((a, (2.0 * norm.sf(np.abs(out[a][1])),
                   bh(2.0 * norm.sf(np.abs(out[a][1]))))) for a in arms)
    hdr = ['domain_idx', 'gene', 'chrom', 'forced', 'n_var_domain',
           'n_var_major', 'n_var_minor', 'n_var_ccre', 'n_individuals',
           'w_major_att']
    with open(OUT + '/att_gene_results.csv', 'w') as fh:
        fh.write(','.join(hdr + sum([['z_' + a, 'p_' + a, 'q_' + a]
                                     for a in arms], [])
                          + ['delta_att_flat', 'delta_int_flat']) + '\n')
        for g in range(C.NDOM):
            row = [g, str(udom[g]), str(dom_chrom[g]), int(dom_forced[g]),
                   int(dn[g]), int(NK[g, 0]), int(NK[g, 1]), int(nccre[g]),
                   int(N), _f(w_att[g])]
            for a in arms:
                row += [_f(out[a][1][g]), _f(qs[a][0][g]), _f(qs[a][1][g])]
            row.append(_f(dz_obs[g]))
            row.append(_f(np.abs(out['INT'][1][g]) - np.abs(zflat[g])))
            fh.write(','.join(str(x) for x in row) + '\n')
    cols = ['log10_nvar', 'size_frac', 'mean_log10_maf', 'mean_r2',
            'frac_re2g', 'frac_ccre']
    with open(OUT + '/att_weights.csv', 'w') as fh:
        fh.write('domain_idx,gene,chrom,w_major_att,w_major_size,'
                 'w_major_att_maf001,'
                 + ','.join('maj_' + x for x in cols) + ','
                 + ','.join('min_' + x for x in cols) + '\n')
        for g in range(C.NDOM):
            fh.write(','.join([str(g), str(udom[g]), str(dom_chrom[g]),
                               _f(w_att[g]), _f(w_size[g]), _f(w_a1[g])]
                              + [_f(RAW[0][g, 0, j]) for j in range(6)]
                              + [_f(RAW[0][g, 1, j]) for j in range(6)]) + '\n')
    with open(OUT + '/att_fit_meta.json', 'w') as fh:
        json.dump(dict(n_individuals=int(N), n_domains=int(C.NDOM),
                       n_domains_modelled=int(gmask.sum()),
                       fold_seed=C.SEED, nfold=NFOLD, ninner=NINNER,
                       lam_grid=LAM, epochs=EPOCHS, lr=LR, hidden=HID,
                       npc=NPC, selection=sel,
                       lam_fixed_for_geneholdout=lam_fix,
                       n_shuffles_run=nshuf, wall_obs_s=round(t_obs, 1),
                       wall_total_s=round(time.time() - t0, 1),
                       covariates=C.COVARS, post_hoc=True, amendment='A5-WIDE',
                       variant=dict(
                           lam_rule='CV-min (no 1-SE relaxation)',
                           gate='MLP forced (hid=%d), w range (%.2f,%.2f)' % (HID, WLO, WHI),
                           base_run='fset/att (1-SE, 2026-09-23 15:02)')),
                  fh, indent=1, default=float, ensure_ascii=False)
    with open(OUT + '/att_fit.done', 'w') as fh:
        fh.write('%.1f s\n' % (time.time() - t0))
    print('DONE %.1f s' % (time.time() - t0))

if __name__ == '__main__':
    main()
