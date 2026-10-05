#!/usr/bin/env python
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import gzip, json, os, sys, bisect, collections, hashlib
import numpy as np
from cyvcf2 import VCF

CHR = sys.argv[1]; GLIST = sys.argv[2]
ROOT = _config_path("${PROJECT_ROOT}/work"); OUT = f"{ROOT}/gate1/out"
ORIG = f"{ROOT}/ref/orig_index/chr{CHR}.vcf.gz"
COMMON = f"{ROOT}/ref/common05/chr{CHR}.maf05.vcf.gz"
CCRE = f"{ROOT}/ref/b6_cards/ccre.s.bed"
SEED = 20260921
VAR_CAP = 1200
NSUB = 2000
RADII_H = [0, 1, 2, 3, 5]
RADII_Q = [0.05, 0.10, 0.25, 0.50]

rng = np.random.default_rng(SEED)

with gzip.open(f"{OUT}/windows_chr{CHR}.json.gz", "rt") as f:
    WIN = json.load(f)
VPOS = {}
with gzip.open(f"{OUT}/vpos_chr{CHR}.tsv.gz", "rt") as f:
    f.readline()
    for line in f:
        a, b, c = line.rstrip("\n").split("\t")
        VPOS[a] = (int(b), int(c))

cc = []
with open(CCRE) as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        if p[0] not in (f"chr{CHR}", CHR): continue
        cc.append((int(p[1]), int(p[2]), p[3] if len(p) > 3 else "."))
cc.sort(); cc_s = [x[0] for x in cc]
def ccre_of(p):
    i = bisect.bisect_right(cc_s, p) - 1
    if i >= 0 and cc[i][1] >= p: return i
    return -1

def parse_id(vid):
    c, pos, ref, alt = vid.split(":")
    return c, int(pos), ref, alt

def read_dosage(ids):
    recs = []
    for vid in ids:
        c, pos, ref, alt = parse_id(vid)
        recs.append((pos, vid, ref, alt))
    recs.sort()
    want = {}
    for pos, vid, ref, alt in recs:
        want.setdefault(pos, []).append((ref, alt, vid))
    regions = []
    for pos in sorted(want):
        if regions and pos - regions[-1][1] < 100000: regions[-1][1] = pos
        else: regions.append([pos, pos])
    out = {}
    for src in (ORIG, COMMON):
        if not os.path.exists(src): continue
        v = VCF(src)
        for s, e in regions:
            try: it = v(f"{CHR}:{max(1,s)}-{e}")
            except Exception: continue
            for rec in it:
                cand = want.get(rec.POS)
                if not cand: continue
                for ref, alt, vid in cand:
                    if vid in out: continue
                    if rec.REF != ref: continue
                    if alt not in (rec.ALT or []): continue
                    ds = rec.format('DS')
                    if ds is None: continue
                    out[vid] = np.asarray(ds, dtype=np.float32).ravel()
        v.close()
    kept = [vid for _, vid, _, _ in recs if vid in out]
    if not kept: return np.zeros((0, 0), dtype=np.float32), []
    X = np.vstack([out[v] for v in kept])
    np.nan_to_num(X, copy=False, nan=0.0)
    return X, kept

def m_eff(Xs):
    if Xs.shape[0] < 2: return float(Xs.shape[0]), 1.0, float('nan')
    sd = Xs.std(axis=1); ok = sd > 0
    if ok.sum() < 2: return float(ok.sum()), 1.0, float('nan')
    Z = (Xs[ok] - Xs[ok].mean(axis=1, keepdims=True)) / sd[ok][:, None]
    R = (Z @ Z.T) / Xs.shape[1]
    ev = np.linalg.eigvalsh(R); ev = np.clip(ev, 0, None)
    M = R.shape[0]
    meff = (ev.sum() ** 2) / (ev ** 2).sum() if (ev ** 2).sum() > 0 else float('nan')
    er = float(np.exp(-(lambda p: (p[p > 0] * np.log(p[p > 0])).sum())(ev / ev.sum())))
    cond = float(ev[-1] / ev[ev > 1e-10][0]) if (ev > 1e-10).any() else float('nan')
    return float(meff), float(meff / M), cond, er, M

def recur_stats(codes):
    _, cnt = np.unique(codes, return_counts=True)
    n = codes.size
    return {
        "unique": int(cnt.size),
        "singleton_frac": float((cnt == 1).sum() / cnt.size),
        "frac_ind_in_ctx_ge2": float(cnt[cnt >= 2].sum() / n),
        "frac_ind_in_ctx_ge5": float(cnt[cnt >= 5].sum() / n),
        "frac_ind_in_ctx_ge10": float(cnt[cnt >= 10].sum() / n),
        "frac_ind_in_ctx_ge50": float(cnt[cnt >= 50].sum() / n),
        "frac_ind_in_ctx_ge100": float(cnt[cnt >= 100].sum() / n),
        "eff_ctx_count": float((cnt.sum() ** 2) / (cnt ** 2).sum()),
        "entropy_nats": float(-(lambda p: (p * np.log(p)).sum())(cnt / cnt.sum())),
        "max_ctx_n": int(cnt.max()),
    }

def hash_rows(B):
    if B.shape[1] == 0: return np.zeros(B.shape[0], dtype=np.int64)
    packed = np.packbits(B, axis=1)
    return np.array([hash(bytes(r)) for r in packed], dtype=np.int64)

def neigh_binary(B, radii):
    if B.shape[1] == 0: return {}
    Bf = B.astype(np.float32)
    s = Bf.sum(axis=1)
    D = s[:, None] + s[None, :] - 2.0 * (Bf @ Bf.T)
    np.fill_diagonal(D, np.inf)
    Ds = np.sort(D, axis=1)
    res = {"nn1": float(np.median(Ds[:, 0])),
           "nn5": float(np.median(Ds[:, 4])) if Ds.shape[1] > 4 else float('nan'),
           "nn10": float(np.median(Ds[:, 9])) if Ds.shape[1] > 9 else float('nan')}
    for r in radii:
        c = (D <= r).sum(axis=1)
        res[f"med_nsim_r{r}_sub"] = float(np.median(c))
        res[f"med_nsim_r{r}_scaled"] = float(np.median(c) * (NSAMP - 1) / (B.shape[0] - 1))
        res[f"frac_ind_nsim_ge50_r{r}"] = float((c * (NSAMP - 1) / (B.shape[0] - 1) >= 50).mean())
    return res

def neigh_cont(E, quants):
    if E.shape[1] == 0: return {}
    Z = (E - E.mean(axis=0)) / (E.std(axis=0) + 1e-9)
    n2 = (Z ** 2).sum(axis=1)
    D = np.maximum(n2[:, None] + n2[None, :] - 2.0 * (Z @ Z.T), 0.0)
    np.fill_diagonal(D, np.inf)
    off = D[np.isfinite(D)]
    res = {}
    Ds = np.sort(D, axis=1)
    res["nn1"] = float(np.median(np.sqrt(Ds[:, 0])))
    res["nn5"] = float(np.median(np.sqrt(Ds[:, 4]))) if Ds.shape[1] > 4 else float('nan')
    for q in quants:
        thr = np.quantile(off, q)
        cval = (D <= thr).sum(axis=1)
        res[f"med_nsim_q{q}_scaled"] = float(np.median(cval) * (NSAMP - 1) / (E.shape[0] - 1))
        res[f"frac_ind_nsim_ge50_q{q}"] = float((cval * (NSAMP - 1) / (E.shape[0] - 1) >= 50).mean())
    return res

genes = [l.strip() for l in open(GLIST) if l.strip()]
NSAMP = None
results = []
for gi, g in enumerate(genes):
    ent = WIN.get(g)
    if ent is None: continue
    bins = {b: [r for r in ent["bins"].get(b, []) if VPOS.get(r[0], (0, 1))[1] == 0] for b in ("B1", "B2", "B3", "B4")}
    nbin_all = {b: len(v) for b, v in bins.items()}
    total = sum(nbin_all.values()); capped = False
    if total > VAR_CAP:
        capped = True
        keep = {}
        for b, v in bins.items():
            k = max(2, int(round(VAR_CAP * len(v) / total))) if len(v) else 0
            k = min(k, len(v))
            idx = rng.choice(len(v), size=k, replace=False) if k < len(v) else np.arange(len(v))
            keep[b] = [v[i] for i in sorted(idx)]
        bins = keep
    ids_by_bin = {b: [r[0] for r in v] for b, v in bins.items()}
    all_ids = [i for b in ("B1", "B2", "B3", "B4") for i in ids_by_bin[b]]
    if len(all_ids) < 2:
        results.append({"gene": g, "status": "too_few_variants", "n_by_bin": nbin_all}); continue
    X, kept = read_dosage(all_ids)
    if X.shape[0] < 2:
        results.append({"gene": g, "status": "dosage_read_failed", "n_by_bin": nbin_all,
                        "n_requested": len(all_ids), "n_read": int(X.shape[0])}); continue
    if NSAMP is None: NSAMP = X.shape[1]
    kept_set = {v: i for i, v in enumerate(kept)}
    binof = {}
    for b in ("B1", "B2", "B3", "B4"):
        for v in ids_by_bin[b]:
            if v in kept_set: binof[kept_set[v]] = b
    B = X >= 0.5
    maf_obs = X.sum(axis=1) / (2.0 * X.shape[1])
    hard_carrier = B.sum(axis=1)
    exp_allele = X.sum(axis=1)
    dvar = X.var(axis=1)
    meta = {r[0]: (r[1], r[2], r[3]) for b in bins for r in bins[b]}
    r2 = np.array([meta[v][1] for v in kept], dtype=np.float32)
    src = np.array([meta[v][2] for v in kept], dtype=np.int8)
    rec = {"gene": g, "chr": CHR, "n_by_bin_full": nbin_all, "capped": capped,
           "n_read": int(X.shape[0]), "nsamp": int(X.shape[1]),
           "window_bp": int(sum(e - s for s, e, _ in ent["span_hg38"])),
           "n_intervals": len(ent["span_hg38"]),
           "src_counts": {"dist_only": int((src == 1).sum()), "re2g_only": int((src == 2).sum()),
                          "both": int((src == 3).sum())}}
    U = {"U1": [i for i in range(X.shape[0]) if binof.get(i) == "B2"],
         "U2": [i for i in range(X.shape[0]) if binof.get(i) in ("B2", "B3")],
         "U3": list(range(X.shape[0]))}
    sub = rng.choice(X.shape[1], size=min(NSUB, X.shape[1]), replace=False)
    per_u = {}
    for uname, idx in U.items():
        if len(idx) == 0: per_u[uname] = {"n_var": 0}; continue
        Xi = X[idx]; Bi = B[idx]
        K = Bi.sum(axis=0); A = Xi.sum(axis=0)
        kd = np.bincount(K, minlength=4)
        carr = int((K >= 1).sum())
        d = {"n_var": len(idx),
             "median_maf": float(np.median(maf_obs[idx])), "iqr_maf": float(np.subtract(*np.percentile(maf_obs[idx], [75, 25]))),
             "median_r2": float(np.median(r2[idx])), "median_dosage_var": float(np.median(dvar[idx])),
             "median_hard_carrier": float(np.median(hard_carrier[idx])),
             "median_expected_allele": float(np.median(exp_allele[idx])),
             "P_K0": float((K == 0).mean()), "P_K1": float((K == 1).mean()),
             "P_Kge2": float((K >= 2).mean()), "P_Kge3": float((K >= 3).mean()),
             "P_Kge2_given_ge1": float((K >= 2).sum() / carr) if carr else float('nan'),
             "P_Kge3_given_ge1": float((K >= 3).sum() / carr) if carr else float('nan'),
             "median_K": float(np.median(K)), "median_A": float(np.median(A))}
        if len(idx) >= 2:
            me = m_eff(Xi)
            d.update({"M_eff": me[0], "M_eff_over_M": me[1], "condition_number": me[2],
                      "eff_rank_entropy": me[3], "M_used": me[4]})
        d["ctxA"] = recur_stats(hash_rows(Bi.T))
        if len(idx) >= 2:
            sd = Xi.std(axis=1); ok = sd > 0
            blocks = np.arange(len(idx))
            if ok.sum() >= 2:
                Z = (Xi[ok] - Xi[ok].mean(axis=1, keepdims=True)) / sd[ok][:, None]
                R2m = ((Z @ Z.T) / Xi.shape[1]) ** 2
                oki = np.where(ok)[0]; lab = {}; nxt = 0
                for a in range(len(oki)):
                    if oki[a] in lab: continue
                    lab[oki[a]] = nxt
                    for bb in range(a + 1, len(oki)):
                        if oki[bb] not in lab and R2m[a, bb] > 0.8: lab[oki[bb]] = nxt
                    nxt += 1
                blocks = np.array([lab.get(i, nxt + i) for i in range(len(idx))])
                nb = len(set(blocks.tolist()))
                Bb = np.zeros((nb, Bi.shape[1]), dtype=bool)
                for bi, lb in enumerate(sorted(set(blocks.tolist()))):
                    Bb[bi] = Bi[blocks == lb].any(axis=0)
                d["n_blocks_r2_0.8"] = int(nb)
                d["ctxB"] = recur_stats(hash_rows(Bb.T))
                d["highLD_pair_frac"] = float((R2m[np.triu_indices_from(R2m, 1)] > 0.8).mean()) if R2m.shape[0] > 1 else float('nan')
        units = collections.defaultdict(list)
        for j in idx:
            p38 = VPOS.get(kept[j], (0, 0))[0]
            u = ccre_of(p38)
            if u >= 0: units[u].append(j)
        d["n_ccre_units"] = len(units)
        d["n_var_in_ccre"] = int(sum(len(v) for v in units.values()))
        if units:
            Ec = np.vstack([X[v].sum(axis=0) for v in units.values()])
            d["ctxC"] = recur_stats(hash_rows((Ec >= 0.5).T))
            Dm = []
            for v in units.values():
                for bb in ("B1", "B2", "B3", "B4"):
                    sel = [j for j in v if binof.get(j) == bb]
                    Dm.append(X[sel].sum(axis=0) if sel else np.zeros(X.shape[1], dtype=np.float32))
            Ed = np.vstack(Dm)
            d["ctxD"] = recur_stats(hash_rows((Ed >= 0.5).T))
            d["nbhd_C"] = neigh_cont(Ec[:, sub].T, RADII_Q)
            d["nbhd_D"] = neigh_cont(Ed[:, sub].T, RADII_Q)
        d["nbhd_A"] = neigh_binary(Bi[:, sub].T, RADII_H)
        if len(idx) >= 2:
            Bf = Bi.astype(np.float32)
            CO = Bf @ Bf.T
            sd = Xi.std(axis=1); ok = sd > 0
            R2p = None
            if ok.sum() >= 2:
                Z = (Xi - Xi.mean(axis=1, keepdims=True)) / (sd[:, None] + 1e-9)
                R2p = ((Z @ Z.T) / Xi.shape[1]) ** 2
            bb = np.array([binof.get(j, "NA") for j in idx])
            pairstats = {}
            for a in range(4):
                for b2 in range(a, 4):
                    ba, bb2 = f"B{a+1}", f"B{b2+1}"
                    ia = np.where(bb == ba)[0]; ib = np.where(bb == bb2)[0]
                    if ia.size == 0 or ib.size == 0: continue
                    sub_co = CO[np.ix_(ia, ib)]
                    if ba == bb2:
                        m = np.triu(np.ones_like(sub_co, dtype=bool), 1)
                        vals = sub_co[m]
                        ld = R2p[np.ix_(ia, ib)][m] if R2p is not None else None
                    else:
                        vals = sub_co.ravel()
                        ld = R2p[np.ix_(ia, ib)].ravel() if R2p is not None else None
                    if vals.size == 0: continue
                    hi = (ld > 0.8) if ld is not None else np.zeros(vals.size, dtype=bool)
                    pairstats[f"{ba}x{bb2}"] = {
                        "n_pairs": int(vals.size), "median": float(np.median(vals)),
                        "iqr": float(np.subtract(*np.percentile(vals, [75, 25]))),
                        "p90": float(np.percentile(vals, 90)), "p95": float(np.percentile(vals, 95)),
                        "p99": float(np.percentile(vals, 99)), "max": float(vals.max()),
                        "n_ge5": int((vals >= 5).sum()), "n_ge10": int((vals >= 10).sum()),
                        "n_ge50": int((vals >= 50).sum()), "n_ge100": int((vals >= 100).sum()),
                        "n_ge500": int((vals >= 500).sum()),
                        "n_ge50_lowLD": int(((vals >= 50) & (~hi)).sum()),
                        "frac_highLD": float(hi.mean())}
            d["cocarrier"] = pairstats
        per_u[uname] = d
    rec["universes"] = per_u
    rec["status"] = "ok"
    results.append(rec)
    if (gi + 1) % 5 == 0:
        print(f"progress chr{CHR} {gi+1}/{len(genes)}", flush=True)

tmp = f"{OUT}/ctx_chr{CHR}.json.tmp"
with open(tmp, "w") as f:
    json.dump({"chr": CHR, "seed": SEED, "var_cap": VAR_CAP, "n_sub_individuals": NSUB,
               "radii_hamming": RADII_H, "radii_quantile": RADII_Q,
               "noncoding_def": "no overlap with gencode protein_coding CDS (hg38)",
               "genes": results}, f)
os.replace(tmp, f"{OUT}/ctx_chr{CHR}.json")
print(json.dumps({"chr": CHR, "genes_done": sum(1 for r in results if r.get("status") == "ok"),
                  "genes_failed": sum(1 for r in results if r.get("status") != "ok")}))
