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
ap = argparse.ArgumentParser(); ap.add_argument("a"); ap.add_argument("basis"); ap.add_argument("--null", type=int, default=100); ap.add_argument("--out", default="")
ap.add_argument("--target-root", default=_config_path("${PROJECT_ROOT}/work/ref15")); ap.add_argument("--chrs", required=True); ap.add_argument("--tag", required=True); A = ap.parse_args()
T0 = time.time(); log = lambda m: print(f"[{time.time()-T0:6.0f}s] {m}", flush=True)
def require(c, m):
    if not c: print("GATE FAIL:", m, flush=True); sys.exit(2)
FOLDS0 = set(A.chrs.split(",")); TRAITS = ["tchl","htn","dm","lip"]
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

R = A.target_root
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
del row_of; NS = N_SAMPLES; log(f"DS loaded {len(DSC)} chr")
sidx = {s: j for j, s in enumerate(SAMPLES)}; TR = {}
for t in TRAITS:
    full = np.full(NS, np.nan)
    with open(_config_path(f'${{PROJECT_ROOT}}/work/ref/annot/offset/{t}.eta.tsv')) as fh:
        hdr = fh.readline().rstrip("\n").split("\t"); iy, ie, im = hdr.index("y"), hdr.index("eta"), hdr.index("mu")
        for line in fh:
            fl = line.rstrip("\n").split("\t"); j = sidx.get(fl[0])
            if j is None: continue
            y = float(fl[iy]); full[j] = y - float(fl[ie]) if t == "tchl" else y - float(fl[im])
    m = ~np.isnan(full); r = full[m].astype(np.float32); TR[t] = dict(cols=np.nonzero(m)[0], r=r - r.mean())
_ord = np.argsort(GID, kind="stable"); _cut = np.searchsorted(GID[_ord], np.arange(NG + 1)); pairs_of = [_ord[_cut[g]:_cut[g+1]] for g in range(NG)]
def burden(w):
    S = np.zeros((NG, NS), np.float32)
    for g in range(NG):
        P = pairs_of[g]; M = DSC[gchr[gl[g]]]; S[g] = (M[VLOC[P]].T @ w[P]).astype(np.float32)
    return S
def zrow(S): return (S - S.mean(1, keepdims=True)) / (S.std(1, keepdims=True) + 1e-8)
def gene_T(S, t):
    c = TR[t]["cols"]; Z = zrow(S[:, c]).astype(np.float64); return (Z @ TR[t]["r"].astype(np.float64)) / (np.sqrt(len(c)) * TR[t]["r"].std())
def mix_terms(T, pi, tau2):
    v = 1.0 + tau2; a_ = tau2 / (2.0 * v); lr = np.exp(np.minimum(T * T * a_, 700.0)) / np.sqrt(v); m = (1 - pi) + pi * lr
    return -np.log(m)
FROZEN = json.load(open(_config_path("${PROJECT_ROOT}/work/ref15/annot/l1/spline_v7_fold0.json")))["frozen_head"]
S_flat = burden(PW); S_learn = burden(phi * PW)
out = dict(tag=A.tag, chrs=sorted(FOLDS0, key=int), genes=NG, pairs=int(len(keys)), clamp_frac=float((np.abs(f) >= CL).mean()), frozen_head=FROZEN, per_trait={})
for t in TRAITS:
    Tf, Tl = gene_T(S_flat, t), gene_T(S_learn, t); pi, tau2 = FROZEN[t]
    lf, ll = mix_terms(Tf, pi, tau2), mix_terms(Tl, pi, tau2); d = lf - ll
    o = np.argsort(-d)
    out["per_trait"][t] = dict(J_flat=float(lf.mean()), J_learn=float(ll.mean()), J_rel=float((ll.mean() - lf.mean()) / abs(lf.mean())),
        ratio=float(np.abs(Tl).mean() / np.abs(Tf).mean()), medT2_flat=float(np.median(Tf**2)), medT2_learn=float(np.median(Tl**2)),
        frac_genes_J_improved=float((d > 0).mean()), top5_share_of_J_gain=float(d[o[:5]].sum() / max(d[d > 0].sum(), 1e-12)), top20_share=float(d[o[:20]].sum() / max(d[d > 0].sum(), 1e-12)),
        n_T_gt3_flat=int((np.abs(Tf) > 3).sum()), n_T_gt3_learn=int((np.abs(Tl) > 3).sum()), maxT_flat=float(np.abs(Tf).max()), maxT_learn=float(np.abs(Tl).max()),
        top5=[dict(gene=gl[i], chr=gchr[gl[i]], T_flat=round(float(Tf[i]), 2), T_learn=round(float(Tl[i]), 2), dJ=round(float(d[i]), 5)) for i in o[:5]],
        worst5=[dict(gene=gl[i], chr=gchr[gl[i]], T_flat=round(float(Tf[i]), 2), T_learn=round(float(Tl[i]), 2), dJ=round(float(d[i]), 5)) for i in o[-5:]])
out["J_flat_mean"] = float(np.mean([out["per_trait"][t]["J_flat"] for t in TRAITS])); out["J_learn_mean"] = float(np.mean([out["per_trait"][t]["J_learn"] for t in TRAITS]))
out["ratio_mean"] = float(np.mean([out["per_trait"][t]["ratio"] for t in TRAITS])); out["sec"] = round(time.time() - T0)
op = A.out or _config_path(f'${{PROJECT_ROOT}}/work/ref15/annot/l1/diag_jratio_{A.tag}.json')
json.dump(out, open(op, "w"), indent=1); print("RESULT " + json.dumps(out)); print("DIAG_DONE")
