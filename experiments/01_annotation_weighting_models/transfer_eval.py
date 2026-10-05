#!/usr/bin/env python3
import os as _cfg_os
import re as _cfg_re
import math as _cfg_math

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)

def _config_number(name, cast=int, positive=False):
    value = cast(_cfg_os.environ[name])
    if not _cfg_math.isfinite(value) or (positive and value <= 0):
        raise ValueError(name + " must be finite" + (" and positive" if positive else ""))
    return value

import os as _os
N_SAMPLES = _config_number("N_SAMPLES", int, True)
import sys, json, time, argparse, pickle, numpy as np, scipy.sparse as sp
from scipy.linalg import null_space
from sklearn.preprocessing import SplineTransformer
ap = argparse.ArgumentParser(); ap.add_argument("a"); ap.add_argument("basis"); ap.add_argument("--null", type=int, default=100); ap.add_argument("--out", default=_config_path("${PROJECT_ROOT}/work/ref15/annot/l1/transfer_fold0.json"))
ap.add_argument("--target-root", default=_config_path("${PROJECT_ROOT}/work/ref")); A = ap.parse_args()
T0 = time.time(); log = lambda m: print(f"[{time.time()-T0:6.0f}s] {m}", flush=True)
def require(c, m):
    if not c: print("GATE FAIL:", m, flush=True); sys.exit(2)
FOLDS0 = {"7","13","16","19","21"}; TRAITS = ["tchl","htn","dm","lip"]
BS = pickle.load(open(A.basis, "rb")); BASIS, CL, coef_names = BS["basis"], BS["CL"], BS["coef_names"]
S1, S2, S3 = BS["cols_S1"], BS["cols_S2"], BS["cols_S3"]
a = np.load(A.a).astype(np.float64); require(len(a) == len(coef_names), f"a len {len(a)} != coef {len(coef_names)}")
R = A.target_root
C = np.load(f"{R}/annot/cache/fm_all.npz", allow_pickle=True)
X, cols, keys, genes, chrs = C["X"], [str(c) for c in C["cols"]], C["key37"].astype(str), C["gene"].astype(str), C["chr"].astype(int).astype(str)
sel = np.nonzero(np.isin(chrs, list(FOLDS0)))[0]
X, keys, genes, chrs = X[sel], keys[sel], genes[sel], chrs[sel]
gl = sorted(set(genes)); gix = {g: i for i, g in enumerate(gl)}; GID = np.array([gix[g] for g in genes]); NG = len(gl)
gchr = {}
for g_, c_ in zip(genes, chrs): gchr.setdefault(g_, c_)
log(f"target band pairs {len(keys):,} genes {NG:,} on FOLDS[0] {sorted(FOLDS0, key=int)}")
def col(cn): return X[:, cols.index(cn)].astype(np.float64)
blocks = []; names = []
def spl_block(cn, x, mask):
    b = BASIS[cn]; xs = np.zeros_like(x); xs[mask] = (x[mask] - b["mu"]) / b["sd"]
    if b["kind"] == "lin":
        blocks.append((xs[:, None] * mask[:, None]).astype(np.float32)); names.append(f"{cn}:0"); return
    kn = np.array(b["knots"]); spl = SplineTransformer(degree=3, knots=kn[:, None], include_bias=True, extrapolation="constant").fit(kn[:, None])
    Braw = spl.transform(xs[:, None]); Q = np.array(b["Q"]); Bt = ((Braw - np.array(b["mean_tr"])) @ Q) * mask[:, None]
    blocks.append(Bt.astype(np.float32)); names.extend(f"{cn}:{i}" for i in range(Bt.shape[1]))
def bin_block(cn, v):
    blocks.append((v - BASIS[cn]["p"])[:, None].astype(np.float32)); names.append(f"{cn}:0")
for cn in S1:
    x = col(cn); mask = ~np.isnan(x)
    if cn == "dist_tss": x = np.log1p(np.nan_to_num(x, nan=0.0)); mask = ~np.isnan(col(cn))
    spl_block(cn, x, mask)
for cn in S2:
    x = col(cn); Amask = ~np.isnan(x)
    if cn == "gh_link_score": x = np.where(Amask, np.log1p(np.nan_to_num(x, nan=0.0)), np.nan)
    spl_block(cn, np.nan_to_num(x, nan=0.0), Amask)
    if cn + "_A" in BASIS: bin_block(cn + "_A", Amask.astype(np.float64))
for cn in S3: bin_block(cn, col(cn))
B = np.hstack(blocks); require(names == coef_names, f"coef order mismatch: {names[:5]} vs {coef_names[:5]}")
f = B.astype(np.float64) @ a; phi = np.exp(np.clip(f, -CL, CL)).astype(np.float32)
PW = (1.0 / np.maximum(col("n_genes"), 1.0)).astype(np.float32)
def oor_frac(cn):
    b = BASIS[cn]; x = col(cn)
    if cn == "dist_tss": x = np.log1p(np.nan_to_num(x, nan=0.0))
    m = ~np.isnan(col(cn)); xs = (x[m] - b["mu"]) / b["sd"]
    return float(((xs < b["knots"][0]) | (xs > b["knots"][-1])).mean())
oor = {cn: oor_frac(cn) for cn in S1 if BASIS[cn]["kind"] == "spl"}
log(f"design {B.shape}; clamp_frac {(np.abs(f) >= CL).mean():.4f}; out-of-knot-range frac per S1 col: {json.dumps({k: round(v,4) for k, v in oor.items()})}")
DSC = {}; VCHR = np.empty(len(keys), dtype=object); VLOC = np.empty(len(keys), np.int64); row_of = {}
for c in sorted(FOLDS0, key=int):
    z = np.load(f"{R}/annot/ds/chr{c}.ds.npz", allow_pickle=True)
    if not DSC: SAMPLES = [str(s_) for s_ in z["samples"]]
    else: require([str(s_) for s_ in z["samples"]] == SAMPLES, f"chr{c} sample order")
    DSC[c] = sp.csr_matrix((z["data"], z["indices"].astype(np.int32, copy=False), z["indptr"].astype(np.int32)), shape=(len(z["keys"]), N_SAMPLES))
    for i_, k in enumerate(z["keys"]): row_of[str(k)] = (c, i_)
    del z
miss = sum(1 for k in keys if str(k) not in row_of); require(miss == 0, f"{miss} pairs lack DS rows")
for p_, k in enumerate(keys): VCHR[p_], VLOC[p_] = row_of[str(k)]
del row_of; NS = N_SAMPLES; log(f"DS loaded {len(DSC)} chr, rows {sum(m.shape[0] for m in DSC.values()):,}")
sidx = {s: j for j, s in enumerate(SAMPLES)}; TR = {}; FULL = {}
for t in TRAITS:
    full = np.full(NS, np.nan)
    with open(f"{R}/annot/offset/{t}.eta.tsv") as fh:
        hdr = fh.readline().rstrip("\n").split("\t"); iy, ie, im = hdr.index("y"), hdr.index("eta"), hdr.index("mu")
        for line in fh:
            fl = line.rstrip("\n").split("\t"); j = sidx.get(fl[0])
            if j is None: continue
            y = float(fl[iy]); full[j] = y - float(fl[ie]) if t == "tchl" else y - float(fl[im])
    FULL[t] = full; m = ~np.isnan(full); r = full[m].astype(np.float32); TR[t] = dict(cols=np.nonzero(m)[0], r=r - r.mean())
_ord = np.argsort(GID, kind="stable"); _cut = np.searchsorted(GID[_ord], np.arange(NG + 1)); pairs_of = [_ord[_cut[g]:_cut[g+1]] for g in range(NG)]
def burden(w):
    S = np.zeros((NG, NS), np.float32)
    for g in range(NG):
        P = pairs_of[g]; M = DSC[gchr[gl[g]]]; S[g] = (M[VLOC[P]].T @ w[P]).astype(np.float32)
    return S
def zrow(S): return (S - S.mean(1, keepdims=True)) / (S.std(1, keepdims=True) + 1e-8)
def gene_T(S, t):
    c = TR[t]["cols"]; Z = zrow(S[:, c]).astype(np.float64); return (Z @ TR[t]["r"].astype(np.float64)) / (np.sqrt(len(c)) * TR[t]["r"].std())
t1 = time.time(); S_flat = burden(PW); S_learn = burden(phi * PW); log(f"burden x2 {time.time()-t1:.0f}s")
ho_l = {t: float(np.abs(gene_T(S_learn, t)).mean()) for t in TRAITS}; ho_f = {t: float(np.abs(gene_T(S_flat, t)).mean()) for t in TRAITS}
ratio = {t: ho_l[t] / ho_f[t] for t in TRAITS}
infl = {t: float(np.median(gene_T(S_learn, t) ** 2) / np.median(gene_T(S_flat, t) ** 2)) for t in TRAITS}
medT2 = {t: [float(np.median(gene_T(S_flat, t) ** 2)), float(np.median(gene_T(S_learn, t) ** 2))] for t in TRAITS}
corr_diag = {cn: float(np.corrcoef(np.nan_to_num(col(cn)), np.log(phi))[0, 1]) for cn in ["maf","r2","avg_cs"]}
null = {}
if A.null:
    rng = np.random.default_rng(123); Zl = {t: zrow(S_learn[:, TR[t]["cols"]]) for t in TRAITS}; Zf = {t: zrow(S_flat[:, TR[t]["cols"]]) for t in TRAITS}
    rat = []
    for b_ in range(A.null):
        perm = rng.permutation(NS); rr = []
        for t in TRAITS:
            full = FULL[t][perm]; m = ~np.isnan(FULL[t]); r_ = full[m]; r_ = np.nan_to_num(r_ - np.nanmean(r_))
            rr.append(np.abs(Zl[t] @ r_).mean() / (np.abs(Zf[t] @ r_).mean() + 1e-12))
        rat.append(rr)
    rat = np.array(rat); obs = np.array([ratio[t] for t in TRAITS]); pm = rat.mean(1)
    null = dict(B=A.null, q95={t: float(np.quantile(rat[:, i], .95)) for i, t in enumerate(TRAITS)}, q95_mean=float(np.quantile(pm, .95)), null_mean={t: float(rat[:, i].mean()) for i, t in enumerate(TRAITS)},
                p={t: float(((rat[:, i] >= obs[i]).sum() + 1) / (A.null + 1)) for i, t in enumerate(TRAITS)}, p_mean=float(((pm >= obs.mean()).sum() + 1) / (A.null + 1)))
out = dict(kind="transfer_1to5pct_phi_applied_to_0.1to1pct", source_a=A.a, basis=A.basis, target_root=R, test_chr=sorted(FOLDS0, key=int), pairs=int(len(keys)), genes=NG,
           clamp_frac=float((np.abs(f) >= CL).mean()), out_of_knot_range=oor, phi_q=[float(q) for q in np.quantile(phi, [0,.01,.5,.99,1])],
           holdout_learned=ho_l, holdout_flat=ho_f, ratio=ratio, ratio_mean=float(np.mean(list(ratio.values()))), inflation_medT2=infl, medT2_flat_learn=medT2,
           corr_logphi_diag=corr_diag, null_fixed_phi=null, sec=round(time.time() - T0))
json.dump(out, open(A.out, "w"), indent=1); print("RESULT " + json.dumps(out)); print("TRANSFER_DONE")
