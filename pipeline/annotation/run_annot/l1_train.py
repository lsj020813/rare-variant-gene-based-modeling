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
import sys, os, json, time, argparse, collections
import numpy as np, scipy.sparse as sp, torch
from scipy.optimize import minimize
from sklearn.preprocessing import SplineTransformer
torch.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", "8")))
R = _config_path("${PROJECT_ROOT}/work/ref")
FOLDS = {0:{"7","13","16","19","21"}, 1:{"2","8","9","10","15"}, 2:{"1","6","14"},
         3:{"3","12","17","22"}, 4:{"4","5","11","18","20"}}
TRAITS = ["tchl","htn","dm","lip"];
LAMS = (None,); LAM_A = 1e-2
BINARY = {"is_cage_prom","in_body","in_tss3kb","in_re2g","t1_na","is_typed","is_indel"}
def require(c, m):
    if not c: raise RuntimeError(f"GATE FAIL: {m}")
ap = argparse.ArgumentParser(); ap.add_argument("--fold", type=int); ap.add_argument("--smoke", type=str)
ap.add_argument("--maxiter", type=int, default=80); ap.add_argument("--gradcheck-only", action="store_true")
ap.add_argument("--gcv-probe", action="store_true"); ap.add_argument("--perm-r", type=int, default=0, help="[v6 누출 대조] >0 이면 시드로 각 형질 잔차 r 을 개인 간 순열 — 유전형-형질 연결 절단"); ap.add_argument("--arm", choices=["spline", "linear", "nn"], default="spline"); ap.add_argument("--out", default=f"{R}/annot/l1")
A = ap.parse_args(); os.makedirs(A.out, exist_ok=True); T0 = time.time()
def log(*x): print(f"[{time.time()-T0:6.0f}s]", *x, flush=True)

C = np.load(f"{R}/annot/cache/fm_all.npz", allow_pickle=True)
X, cols, keys, genes, chrs = C["X"], [str(c) for c in C["cols"]], C["key37"].astype(str), C["gene"].astype(str), C["chr"].astype(int).astype(str)
require(X.shape == (_config_number("N_ANNOTATION_VARIANTS", int, True), 33), f"cache shape {X.shape}")
if A.smoke:
    sel = np.nonzero(chrs == A.smoke)[0]; use_chr = {A.smoke}
else:
    sel = np.arange(len(keys)); use_chr = set(str(c) for c in range(1, 23))
X, keys, genes, chrs = X[sel], keys[sel], genes[sel], chrs[sel]
gl = sorted(set(genes)); gix = {g: i for i, g in enumerate(gl)}; GID = np.array([gix[g] for g in genes]); NG = len(gl)
gchr = {g: c for g, c in zip(genes, chrs)}
log(f"pairs {len(keys):,} genes {NG:,} chr {sorted(use_chr, key=int)} device {torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'cpu'}")

metas = []
for c in sorted(use_chr, key=int):
    z = np.load(f"{R}/annot/ds/chr{c}.ds.npz", allow_pickle=True, mmap_mode=None)
    require(len(z["samples"]) == N_SAMPLES, "samples"); metas.append((c, len(z["keys"]), int(z["indptr"][-1])))
    if len(metas) == 1: SAMPLES = [str(s_) for s_ in z["samples"]]
    del z
NV, NNZ = sum(m[1] for m in metas), sum(m[2] for m in metas)
data = np.empty(NNZ, np.float32); indices = np.empty(NNZ, np.int32); indptr = np.empty(NV + 1, np.int64); indptr[0] = 0
row_of = {}; off = 0; noff = 0
for c, nv, nnz in metas:
    z = np.load(f"{R}/annot/ds/chr{c}.ds.npz", allow_pickle=True)
    data[noff:noff+nnz] = z["data"]; indices[noff:noff+nnz] = z["indices"]; indptr[off+1:off+nv+1] = z["indptr"][1:] + noff
    for i, k in enumerate(z["keys"]): row_of[str(k)] = off + i
    off += nv; noff += nnz; del z
DS = sp.csr_matrix((data, indices, indptr), shape=(NV, N_SAMPLES)); del data, indices, indptr
VROW = np.array([row_of[str(k)] for k in keys]); require(len(VROW) == len(keys), "pair→DS row")
NS = N_SAMPLES; log(f"DS {DS.shape} nnz {DS.nnz:,} ({DS.nnz/DS.shape[0]/NS*100:.2f}%)")

sidx = {s: j for j, s in enumerate(SAMPLES)}; TR = {}
for t in TRAITS:
    cols_t, r = [], []
    with open(f"{R}/annot/offset/{t}.eta.tsv") as fh:
        hdr = fh.readline().rstrip("\n").split("\t"); iy, ie, im = hdr.index("y"), hdr.index("eta"), hdr.index("mu")
        for line in fh:
            f = line.rstrip("\n").split("\t"); j = sidx.get(f[0])
            if j is None: continue
            y = float(f[iy]); res = y - float(f[ie]) if t == "tchl" else y - float(f[im])
            cols_t.append(j); r.append(res)
    r = np.array(r, np.float32)
    if A.perm_r: r = np.random.default_rng(A.perm_r + len(TR)).permutation(r); log(f"  [PERM-R] {t}: residual permuted across individuals (seed {A.perm_r})")
    TR[t] = dict(cols=np.array(cols_t), r=r - r.mean(), var=float(r.var()))
    require(len(cols_t) > 70000, f"{t} n={len(cols_t)}"); log(f"{t}: n={len(cols_t):,} var={TR[t]['var']:.4f}")

_ord = np.argsort(GID, kind="stable"); _cut = np.searchsorted(GID[_ord], np.arange(NG + 1))
pairs_of = [_ord[_cut[g]:_cut[g+1]] for g in range(NG)]
require(sum(len(p) for p in pairs_of) == len(GID), "pairs_of")
def burden(phi):
    S = np.zeros((NG, NS), phi.dtype)
    for g in range(NG):
        P = pairs_of[g]; S[g] = (DS[VROW[P]].T @ phi[P]).astype(phi.dtype)
    return S
def grad_phi(gS):
    gm = np.zeros(len(VROW), gS.dtype)
    for g in range(NG):
        P = pairs_of[g]; gm[P] = DS[VROW[P]] @ gS[g]
    return gm
def zrow(S): return (S - S.mean(1, keepdims=True)) / (S.std(1, keepdims=True) + 1e-8)
DEV = "cuda" if torch.cuda.is_available() else "cpu"
DT = [torch.float32]
LAM_USED = {}
GCV_LOGMAX = float(os.environ.get("GCV_LOGMAX", "9"))
GCV_GRID = [float(x) for x in np.logspace(-2, GCV_LOGMAX, int(4 * (GCV_LOGMAX + 2)) + 1)]
GCV_TRACE = os.environ.get("GCV_TRACE") == "1"
def gcv_lambda(G, Zr, r, n):
    d, Q = torch.linalg.eigh(G.float()); d = d.double().clamp_min(0.0)
    w = (Q.T @ Zr.float()).double(); rr = float((r.double() @ r.double()).item()); del Q
    if DEV == "cuda": torch.cuda.empty_cache()
    best = (float("inf"), None)
    for lam in GCV_GRID:
        den = d + lam; trH = float((d / den).sum().item())
        rss = rr - float((w * w * (2.0 / den - d / den ** 2)).sum().item())
        g = (rss / n) / max(1e-12, (1.0 - trH / n)) ** 2
        if GCV_TRACE: log(f"    gcv λ={lam:.3g} trH={trH:.1f} rss/n={rss/n:.6f} gcv={g:.6f}")
        if g < best[0]: best = (g, lam)
    return best[1]

HEAD_MODE = os.environ.get("HEAD", "mix")
def mix_em(T, iters=200):
    T2 = T ** 2; pi, v = 0.05, 4.0
    for _ in range(iters):
        f0 = np.exp(-T2 / 2); f1 = np.exp(-T2 / (2 * v)) / np.sqrt(v)
        w = pi * f1 / ((1 - pi) * f0 + pi * f1 + 1e-300)
        pi_n = float(np.clip(w.mean(), 1e-3, 0.5)); v_n = float(max(1.01, (w * T2).sum() / max(w.sum(), 1e-12)))
        if abs(pi_n - pi) < 1e-7 and abs(v_n - v) < 1e-6: pi, v = pi_n, v_n; break
        pi, v = pi_n, v_n
    return pi, v - 1.0

def mix_obj(T, pi, tau2):
    v = 1.0 + tau2; a = tau2 / (2.0 * v)
    lr = np.exp(np.minimum(T * T * a, 700.0)) / np.sqrt(v)
    m = (1 - pi) + pi * lr; G = len(T)
    J = -np.log(m).sum() / G
    w = pi * lr / m
    dJdT = -(w * T * 2 * a) / G
    return J, dJdT, w

def head(S, train_gi, lam, want_grad=True, betas=None):
    npdt = np.float64 if DT[0] == torch.float64 else np.float32
    gS = np.zeros(S.shape, npdt) if want_grad else None; J_tot = 0.0; out = {}
    for t in TRAITS:
        c = TR[t]["cols"]; n = len(c)
        Z = torch.tensor(S[np.ix_(train_gi, c)], device=DEV, dtype=DT[0])
        r = torch.tensor(TR[t]["r"], device=DEV, dtype=DT[0])
        mu = Z.mean(1, keepdim=True); sd = Z.std(1, keepdim=True) + 1e-8
        Z.sub_(mu).div_(sd)
        if HEAD_MODE == "mix":
            sig = float(np.std(TR[t]["r"])); r64 = r.double(); Tt = torch.empty(Z.shape[0], dtype=torch.float64, device=DEV)
            for b0 in range(0, Z.shape[0], 2048):
                Tt[b0:b0 + 2048] = Z[b0:b0 + 2048].double() @ r64
            T = (Tt / (np.sqrt(n) * sig)).cpu().numpy(); del Tt, r64
            if betas is None: pi, tau2 = mix_em(T)
            else: pi, tau2 = betas[t]
            J, dJdT, w = mix_obj(T, pi, tau2)
            out[t] = (pi, tau2); LAM_USED[t] = (round(pi, 5), round(tau2, 4), int((w > 0.5).sum())); J_tot += float(J)
            if want_grad:
                coef = torch.tensor(dJdT / (np.sqrt(n) * sig), device=DEV, dtype=DT[0])
                gZ = coef[:, None] * r[None, :]
                gZ.sub_(gZ.mean(1, keepdim=True)).sub_(Z * (gZ * Z).mean(1, keepdim=True)).div_(sd)
                gS[np.ix_(train_gi, c)] += gZ.cpu().numpy(); del gZ, coef
            del Z, r
            if DEV == "cuda": torch.cuda.empty_cache()
            continue
        if betas is None:
            G = Z @ Z.T; Zr = Z @ r
            if lam is None:
                lam_t = gcv_lambda(G, Zr, r, n)
            else:
                lam_t = float(lam)
            G.diagonal().add_(lam_t); beta = torch.linalg.solve(G, Zr); del G, Zr
        else:
            beta, lam_t = betas[t]; beta = beta.to(DEV, DT[0])
        out[t] = (beta.detach().to("cpu"), lam_t); LAM_USED[t] = lam_t
        e = r - beta @ Z
        e64 = e.double(); b64 = beta.double()
        J = (e64 @ e64 + lam_t * (b64 @ b64)) / (n * TR[t]["var"]); J_tot += float(J.item()); del e64, b64
        if want_grad:
            gZ = (-2.0 / n / TR[t]["var"]) * beta[:, None] * e[None, :]
            gZ.sub_(gZ.mean(1, keepdim=True)).sub_(Z * (gZ * Z).mean(1, keepdim=True)).div_(sd)
            gS[np.ix_(train_gi, c)] += gZ.cpu().numpy(); del gZ
        del Z, e, beta, r
        if DEV == "cuda": torch.cuda.empty_cache()
    return J_tot, gS, out

def holdout(S, test_gi):
    out = {}
    for t in TRAITS:
        c, r = TR[t]["cols"], TR[t]["r"]; Z = zrow(S[test_gi][:, c])
        out[t] = float(np.abs(Z @ r).mean() / (len(r) * r.std() + 1e-8))
    return out

if A.smoke:
    rng = np.random.default_rng(0); perm = rng.permutation(NG); test_gi = np.sort(perm[:NG//5]); train_gi = np.sort(perm[NG//5:])
else:
    test_gi = np.array([gix[g] for g in gl if gchr[g] in FOLDS[A.fold]]); train_gi = np.array([gix[g] for g in gl if gchr[g] not in FOLDS[A.fold]])
require(len(test_gi) and len(train_gi) and not set(test_gi) & set(train_gi), "fold split")
train_pairs = np.nonzero(np.isin(GID, train_gi))[0]
log(f"train genes {len(train_gi):,} test {len(test_gi):,}")

cont = [j for j, cname in enumerate(cols) if cname not in BINARY]; binc = [j for j, cname in enumerate(cols) if cname in BINARY]
Xf = X.copy(); med = np.nanmedian(X[train_pairs], axis=0); nanm = np.isnan(Xf)
Xf[nanm] = np.take(med, np.nonzero(nanm)[1])
mu, sd = Xf[train_pairs].mean(0), Xf[train_pairs].std(0) + 1e-8
Xs = (Xf - mu) / sd
if A.arm == "spline":
    spl = SplineTransformer(n_knots=6, degree=3, knots="quantile", include_bias=False).fit(Xs[train_pairs][:, cont])
    B = np.hstack([spl.transform(Xs[:, cont]), Xs[:, binc]]).astype(np.float32)
    log(f"basis {B.shape}  (cont {len(cont)}×{spl.n_features_out_//len(cont)} + bin {len(binc)})")
else:
    B = Xs.astype(np.float32)
    log(f"basis {B.shape}  ({A.arm})")
if A.arm == "nn":
    torch.manual_seed(0)
    NET = torch.nn.Sequential(torch.nn.Linear(B.shape[1], 64), torch.nn.GELU(), torch.nn.Linear(64, 32), torch.nn.GELU(), torch.nn.Linear(32, 1)).to(DEV)
    NPAR = sum(p.numel() for p in NET.parameters()); log(f"nn params {NPAR}")
    Bt = torch.tensor(B, device=DEV)
    def nn_pack():   return torch.cat([p.detach().reshape(-1) for p in NET.parameters()]).cpu().numpy().astype(np.float64)
    def nn_unpack(v):
        i = 0
        for p in NET.parameters():
            n_ = p.numel(); p.data.copy_(torch.tensor(v[i:i+n_], dtype=p.dtype, device=DEV).reshape(p.shape)); i += n_
    def nn_phi(): return torch.exp(NET(Bt).squeeze(1).clamp(max=6.0))

NIT = [0]
BETAS = [None]
def objective(a, lam):
    dt = np.float64 if DT[0] == torch.float64 else np.float32
    if A.arm == "nn":
        nn_unpack(a); NET.zero_grad()
        with torch.no_grad(): phi = nn_phi().cpu().numpy().astype(dt)
        S = burden(phi).astype(dt); loss, gS, _ = head(S, train_gi, lam, betas=BETAS[0])
        gm_np = grad_phi(gS)
        phi_t = nn_phi(); gm = torch.tensor(gm_np, device=DEV, dtype=phi_t.dtype)
        phi_t.backward(gradient=gm); del phi_t, gm
        ga = torch.cat([p.grad.reshape(-1) for p in NET.parameters()]).cpu().numpy().astype(np.float64) + LAM_A * a * 2 / len(a)
    else:
        a = a.astype(dt); eta = B.astype(dt) @ a
        inside = (np.abs(eta) < 6.0); phi = np.exp(np.clip(eta, -6.0, 6.0)).astype(dt)
        S = burden(phi).astype(dt); loss, gS, _ = head(S, train_gi, lam, betas=BETAS[0])
        if not np.isfinite(loss):
            log(f"  non-finite loss at |a|={np.abs(a).max():.2f} — returning 1e6"); return 1e6, (LAM_A * a * 2 / len(a)).astype(np.float64)
        gm = grad_phi(gS) * phi * inside
        ga = B.T @ gm + LAM_A * a * 2 / len(a)
    NIT[0] += 1
    if NIT[0] % 5 == 1: log(f"  it {NIT[0]} loss {loss:.5f} |a|={np.abs(a).max():.3f}")
    return float(loss + LAM_A * (a ** 2).mean()), np.asarray(ga, np.float64)
def phi_of(a):
    if A.arm == "nn": nn_unpack(a); return nn_phi().detach().cpu().numpy().astype(np.float32)
    return np.exp(np.clip(B @ a.astype(np.float32), -6.0, 6.0)).astype(np.float32)
def fit(lam, maxiter, outer_max=8, inner=15, rel_tol=1e-5):
    a = (nn_pack() if A.arm == "nn" else np.zeros(B.shape[1])).astype(np.float64); NIT[0] = 0
    prev = None; res = None; start = 0
    ck = f"{A.out}/{CKTAG}.ckpt.npz"
    if os.path.exists(ck):
        z = np.load(ck, allow_pickle=True)
        if int(z["n_train"]) == len(train_gi) and z["a"].shape == a.shape:
            a = z["a"].astype(np.float64); start = int(z["outer"]) + 1; prev = float(z["prev"]) if np.isfinite(z["prev"]) else None; NIT[0] = int(z["nit"])
            log(f"  RESUME from {ck}: outer {start} |a|max={np.abs(a).max():.3f} nit {NIT[0]}")
        else: log(f"  ckpt ignored (shape/train mismatch)")
    for outer in range(start, outer_max):
        BETAS[0] = None
        S = burden(phi_of(a)); J0, _, betas = head(S, train_gi, lam, want_grad=False); BETAS[0] = betas
        log(f"  outer {outer} J(β̂ 재적합) {J0:.6f} |a|max={np.abs(a).max():.3f}")
        if prev is not None and (prev - J0) < rel_tol * abs(prev):
            log(f"  outer stop: rel gain {(prev - J0)/abs(prev):.2e} < {rel_tol}"); break
        prev = J0
        res = minimize(objective, a, args=(lam,), jac=True, method="L-BFGS-B", options=dict(maxiter=inner, maxfun=3 * inner, ftol=1e-8, gtol=1e-10, maxcor=20))
        a = res.x.astype(np.float64)
        np.savez(ck + ".tmp.npz", a=a, outer=outer, prev=prev, nit=NIT[0], n_train=len(train_gi)); os.replace(ck + ".tmp.npz", ck)
        if NIT[0] >= maxiter: log(f"  maxiter {maxiter} reached"); break
    BETAS[0] = None
    return a.astype(np.float32), res

rng = np.random.default_rng(1 + (A.fold or 0)); p2 = rng.permutation(train_gi); inner_va, inner_tr = np.sort(p2[:len(p2)//5]), np.sort(p2[len(p2)//5:])
t1 = time.time(); S_flat = burden(np.ones(len(VROW), np.float32)); log(f"burden(flat) {time.time()-t1:.1f}s")
t1 = time.time(); _l, _g, _b = head(S_flat, train_gi, None); log(f"head[{HEAD_MODE}] params (flat): {LAM_USED}"); log(f"head {time.time()-t1:.1f}s"); t1 = time.time(); _ = grad_phi(_g); log(f"grad_phi {time.time()-t1:.1f}s")
if A.gcv_probe: print("GCV_PROBE_DONE"); sys.exit(0)
if A.smoke:
    rng0 = np.random.default_rng(7); npar = NPAR if A.arm == "nn" else B.shape[1]
    a_t = (nn_pack() if A.arm == "nn" else np.zeros(npar)) + rng0.normal(0, 0.05, npar); d = rng0.normal(0, 1, npar); d /= np.linalg.norm(d)
    DT[0] = torch.float64
    if A.arm == "nn": NET.double(); Bt = Bt.double()
    _S = burden(phi_of(a_t.astype(np.float32)).astype(np.float64)); _, _, BETAS[0] = head(_S, train_gi, None, want_grad=False)
    f0, g0 = objective(a_t, None); eps = 1e-2; f1, _ = objective(a_t + eps * d, 10.0); f2, _ = objective(a_t - eps * d, 10.0)
    DT[0] = torch.float32; BETAS[0] = None
    if A.arm == "nn": NET.float(); Bt = Bt.float()
    fd, an = (f1 - f2) / (2 * eps), float(g0 @ d); log(f"GRADCHECK fd={fd:.6e} analytic={an:.6e} rel={abs(fd-an)/(abs(fd)+1e-12):.3e}")
    require(abs(fd - an) / (abs(fd) + 1e-12) < 0.05, "gradient check failed")
    NIT[0] = 0
    if A.gradcheck_only: print("GRADCHECK_ONLY_DONE"); sys.exit(0)
CKTAG = A.arm + "_" + ("smoke" if A.smoke else f"fold{A.fold}") + (f"_chr{A.smoke}" if A.smoke else "") + (f"_perm{A.perm_r}" if A.perm_r else "")
sel_scores = {}
LAM = None
a, res = fit(LAM, maxiter=A.maxiter); log(f"head[{HEAD_MODE}] params (final outer): {LAM_USED}")
S_learn = burden(phi_of(a))
ho_l, ho_f = holdout(S_learn, test_gi), holdout(S_flat, test_gi)
ratio = {t: ho_l[t] / ho_f[t] for t in TRAITS}
out = dict(arm=A.arm, mode="smoke" if A.smoke else f"fold{A.fold}", chr=sorted(use_chr, key=int), pairs=int(len(keys)), genes=NG,
           train_genes=int(len(train_gi)), test_genes=int(len(test_gi)), basis=int(B.shape[1]), head=HEAD_MODE, perm_r=A.perm_r, head_params={t: (list(v) if isinstance(v, tuple) else float(v)) for t, v in LAM_USED.items()}, lam_a=LAM_A, lam_a_form="mean",
           nit=int(res.nit), converged=bool(res.success), loss=float(res.fun), holdout_learned=ho_l, holdout_flat=ho_f, ratio=ratio,
           ratio_mean=float(np.mean(list(ratio.values()))), sec=round(time.time() - T0), inner=sel_scores)
tag = CKTAG
np.save(f"{A.out}/{tag}.a.npy", a); np.save(f"{A.out}/{tag}.phi.npy", phi_of(a)); json.dump(out, open(f"{A.out}/{tag}.json", "w"), indent=1)
ck = f"{A.out}/{CKTAG}.ckpt.npz"; os.path.exists(ck) and os.remove(ck)
log("RESULT", json.dumps(out)); print("L1_DONE")
