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
import sys, os, json, time, argparse, hashlib
import numpy as np, scipy.sparse as sp, torch
from scipy.optimize import minimize
from scipy.linalg import null_space
from sklearn.preprocessing import SplineTransformer
torch.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", "8")))
R = _config_path("${PROJECT_ROOT}/work/ref")
CODE_SHA = hashlib.sha256(open(__file__, "rb").read()).hexdigest()
FOLDS = {0:{"7","13","16","19","21"}, 1:{"2","8","9","10","15"}, 2:{"1","6","14"},
         3:{"3","12","17","22"}, 4:{"4","5","11","18","20"}}
TRAITS = ["tchl","htn","dm","lip"]
S1 = ["cons","epi_active","epi_repr","epi_trans","tf","linsight","gpn_msa","dist_tss"]
S2 = ["re2g_max","gh_elem_score","gh_link_score"]
S3 = ["is_cage_prom","in_body","in_tss3kb","in_re2g","t1_na","is_indel"]
DIAG = ["maf","r2","avg_cs"]
CL = 2.0
EPS_F = 0.05
TIE = 1e-3
NKNOT = 6
def require(c, m):
    if not c: raise RuntimeError(f"GATE FAIL: {m}")
ap = argparse.ArgumentParser(); ap.add_argument("--fold", type=int); ap.add_argument("--smoke", type=str)
ap.add_argument("--arm", default="spline", choices=["spline","linear","nn"]); ap.add_argument("--out", default=f"{R}/annot/l1")
ap.add_argument("--maxiter", type=int, default=60); ap.add_argument("--inner-maxiter", type=int, default=30)
ap.add_argument("--gradcheck-only", action="store_true"); ap.add_argument("--perm-r", type=int, default=0, help="개인 공동 순열 시드 (4형질 동일)")
ap.add_argument("--null-perm", type=int, default=0, help="고정 φ 순열 귀무 횟수 (시험 유전자, 학습 없음)")
ap.add_argument("--lam", type=float, default=None, help="λ 고정 (격자 생략)"); ap.add_argument("--lam-grid", type=int, default=3)
A = ap.parse_args(); os.makedirs(A.out, exist_ok=True); T0 = time.time()
def log(*x): print(f"[{time.time()-T0:6.0f}s]", *x, flush=True)

C = np.load(f"{R}/annot/cache/fm_all.npz", allow_pickle=True)
X, cols, keys, genes, chrs = C["X"], [str(c) for c in C["cols"]], C["key37"].astype(str), C["gene"].astype(str), C["chr"].astype(int).astype(str)
require(X.shape == (_config_number("N_ANNOTATION_VARIANTS", int, True), 33), f"cache shape {X.shape}")
for cn in S1 + S2 + S3 + DIAG: require(cn in cols, f"column {cn} missing from cache")
if A.smoke:
    sel = np.nonzero(chrs == A.smoke)[0]; use_chr = {A.smoke}
else:
    sel = np.arange(len(keys)); use_chr = set(str(c) for c in range(1, 23))
X, keys, genes, chrs = X[sel], keys[sel], genes[sel], chrs[sel]
gl = sorted(set(genes)); gix = {g: i for i, g in enumerate(gl)}; GID = np.array([gix[g] for g in genes]); NG = len(gl)
gchr = {g: c for g, c in zip(genes, chrs)}
log(f"pairs {len(keys):,} genes {NG:,} chr {sorted(use_chr, key=int)} device {torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'cpu'} code {CODE_SHA[:12]}")

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
FULL = {}
for t in TRAITS:
    full = np.full(NS, np.nan, np.float64)
    with open(f"{R}/annot/offset/{t}.eta.tsv") as fh:
        hdr = fh.readline().rstrip("\n").split("\t"); iy, ie, im = hdr.index("y"), hdr.index("eta"), hdr.index("mu")
        for line in fh:
            f = line.rstrip("\n").split("\t"); j = sidx.get(f[0])
            if j is None: continue
            y = float(f[iy]); full[j] = y - float(f[ie]) if t == "tchl" else y - float(f[im])
    FULL[t] = full
if A.perm_r:
    perm = np.random.default_rng(A.perm_r).permutation(NS)
    for t in TRAITS: FULL[t] = FULL[t][perm]
    log(f"  [PERM-R] joint individual permutation seed {A.perm_r}")
for t in TRAITS:
    m = ~np.isnan(FULL[t]); r = FULL[t][m].astype(np.float32)
    TR[t] = dict(cols=np.nonzero(m)[0], r=r - r.mean(), var=float(r.var()))
    require(len(r) > 70000, f"{t} n={len(r)}"); log(f"{t}: n={len(r):,} var={TR[t]['var']:.4f}")

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

def head(S, gene_idx, frozen, want_grad=True):
    npdt = np.float64 if DT[0] == torch.float64 else np.float32
    gS = np.zeros(S.shape, npdt) if want_grad else None; J_tot = 0.0
    for t in TRAITS:
        c = TR[t]["cols"]; n = len(c)
        Z = torch.tensor(S[np.ix_(gene_idx, c)], device=DEV, dtype=DT[0]); r = torch.tensor(TR[t]["r"], device=DEV, dtype=DT[0])
        mu = Z.mean(1, keepdim=True); sd = Z.std(1, keepdim=True) + 1e-8; Z.sub_(mu).div_(sd)
        sig = float(np.std(TR[t]["r"])); r64 = r.double(); Tt = torch.empty(Z.shape[0], dtype=torch.float64, device=DEV)
        for b0 in range(0, Z.shape[0], 2048): Tt[b0:b0 + 2048] = Z[b0:b0 + 2048].double() @ r64
        T = (Tt / (np.sqrt(n) * sig)).cpu().numpy(); del Tt, r64
        pi, tau2 = frozen[t]; J, dJdT, w = mix_obj(T, pi, tau2); J_tot += float(J)
        if want_grad:
            coef = torch.tensor(dJdT / (np.sqrt(n) * sig), device=DEV, dtype=DT[0]); gZ = coef[:, None] * r[None, :]
            gZ.sub_(gZ.mean(1, keepdim=True)).sub_(Z * (gZ * Z).mean(1, keepdim=True)).div_(sd)
            gS[np.ix_(gene_idx, c)] += gZ.cpu().numpy(); del gZ, coef
        del Z, r
        if DEV == "cuda": torch.cuda.empty_cache()
    return J_tot, gS
def freeze_head(S, gene_idx):
    fr = {}
    for t in TRAITS:
        c = TR[t]["cols"]; Z = zrow(S[np.ix_(gene_idx, c)]).astype(np.float64); T = (Z @ TR[t]["r"].astype(np.float64)) / (np.sqrt(len(c)) * TR[t]["r"].std())
        fr[t] = mix_em(T)
    return fr
def gene_T(S, gene_idx, t):
    c = TR[t]["cols"]; Z = zrow(S[np.ix_(gene_idx, c)]).astype(np.float64)
    return (Z @ TR[t]["r"].astype(np.float64)) / (np.sqrt(len(c)) * TR[t]["r"].std())
def holdout(S, test_gi):
    return {t: float(np.abs(gene_T(S, test_gi, t)).mean()) for t in TRAITS}

if A.smoke:
    rng = np.random.default_rng(0); perm = rng.permutation(NG); test_gi = np.sort(perm[:NG//5]); train_gi = np.sort(perm[NG//5:])
    p2 = np.random.default_rng(1).permutation(train_gi); inner_va, inner_tr = np.sort(p2[:len(p2)//5]), np.sort(p2[len(p2)//5:])
else:
    test_gi = np.array([gix[g] for g in gl if gchr[g] in FOLDS[A.fold]]); train_gi = np.array([gix[g] for g in gl if gchr[g] not in FOLDS[A.fold]])
    VA_CHR = FOLDS[(A.fold + 1) % 5]
    inner_va = np.array([g for g in train_gi if gchr[gl[g]] in VA_CHR]); inner_tr = np.array([g for g in train_gi if gchr[gl[g]] not in VA_CHR])
require(len(test_gi) and len(train_gi) and not set(test_gi) & set(train_gi), "fold split")
require(len(inner_va) and len(inner_tr) and not set(inner_va) & set(inner_tr), "inner split")
train_pairs = np.nonzero(np.isin(GID, train_gi))[0]
log(f"train genes {len(train_gi):,} (inner tr {len(inner_tr):,} / va {len(inner_va):,}) test {len(test_gi):,}")

def col(cn): return X[:, cols.index(cn)].astype(np.float64)
BLOCKS = []
IS_TRAIN = np.zeros(len(keys), bool); IS_TRAIN[train_pairs] = True
def std_obs(x, mask):
    tr = mask & IS_TRAIN
    mu, sd = x[tr].mean(), x[tr].std() + 1e-8
    xs = np.zeros_like(x); xs[mask] = (x[mask] - mu) / sd
    return xs, tr
def add_linear_cont(cn, x, mask):
    xs, tr = std_obs(x, mask)
    BLOCKS.append((cn, "lin", xs[:, None].astype(np.float32), np.array([[float((xs[train_pairs] ** 2).mean())]])))
def add_spline_cont(cn, x, mask):
    xs, tr = std_obs(x, mask)
    qk = np.quantile(xs[tr], np.linspace(0, 1, NKNOT)); kn = [qk[0]]
    for v in qk[1:]:
        if v - kn[-1] > 1e-6 * (qk[-1] - qk[0] + 1e-12): kn.append(v)
    kn = np.array(kn)
    if len(kn) < 3:
        add_linear_cont(cn, x, mask); return dict(fallback="linear", unique_knots=int(len(kn)))
    spl = SplineTransformer(degree=3, knots=kn[:, None], include_bias=True, extrapolation="constant").fit(xs[tr][:, None])
    Braw = spl.transform(xs[:, None])
    mean_tr = Braw[tr].mean(0)
    Q = null_space(np.ones((1, Braw.shape[1])))
    Bt = ((Braw - mean_tr) @ Q) * mask[:, None]
    G = Bt[tr].T @ Bt[tr] / tr.sum()
    bs = spl.bsplines_[0]; tt = bs.t; lo, hi = tt[3], tt[-4]; grid = np.linspace(lo, hi, 2001)
    D2 = bs.derivative(2)(grid)
    Om_raw = D2.T @ D2 * (grid[1] - grid[0]); Om = Q.T @ Om_raw @ Q
    ev = np.linalg.eigvalsh(Om); rank = int((ev > 1e-8 * ev.max()).sum())
    q = float(np.trace(np.linalg.solve(G, Om)) / rank)
    P = G + Om / q
    BLOCKS.append((cn, "spl", Bt.astype(np.float32), P))
    return dict(rank_omega=rank, q=q, n_basis=int(Braw.shape[1]), knots=[round(float(k), 3) for k in kn], n_obs_train=int(tr.sum()))
def add_binary(cn, b, kind="bin"):
    p = float(b[train_pairs].mean())
    BLOCKS.append((cn, kind, (b - p)[:, None].astype(np.float32), np.array([[1.0]])))
    return p
META = {}
for cn in S1:
    x = col(cn); mask = ~np.isnan(x)
    if cn == "dist_tss": x = np.log1p(np.nan_to_num(x, nan=0.0)); mask = ~np.isnan(col(cn))
    META[cn] = (add_spline_cont(cn, x, mask) if A.arm == "spline" else add_linear_cont(cn, x, mask)) or {}
    META[cn]["na_pct"] = round(100 * (1 - mask.mean()), 2)
BIN = {cn: col(cn) for cn in S3}
for cn in S2:
    x = col(cn); Amask = ~np.isnan(x)
    if cn == "gh_link_score": x = np.where(Amask, np.log1p(np.nan_to_num(x, nan=0.0)), np.nan)
    META[cn] = (add_spline_cont(cn, np.nan_to_num(x, nan=0.0), Amask) if A.arm == "spline" else add_linear_cont(cn, np.nan_to_num(x, nan=0.0), Amask)) or {}
    META[cn]["applicable_pct"] = round(100 * Amask.mean(), 2)
    dup = [b for b in S3 if np.array_equal(BIN[b].astype(bool), Amask)]
    if dup: META[cn]["indicator"] = f"dup of {dup[0]} — skipped"
    else:   META[cn]["indicator"] = f"delta p={add_binary(cn + '_A', Amask.astype(np.float64), 'ind'):.4f}"
for cn in S3: META[cn] = dict(p=add_binary(cn, BIN[cn]))
B = np.hstack([b[2] for b in BLOCKS]); COEF_NAMES = [f"{b[0]}:{i}" for b in BLOCKS for i in range(b[2].shape[1])]
Pfull = np.zeros((B.shape[1], B.shape[1])); i0 = 0
for b in BLOCKS:
    d = b[2].shape[1]; Pfull[i0:i0+d, i0:i0+d] = b[3]; i0 += d
log(f"design {B.shape} blocks {len(BLOCKS)} ({A.arm}); S2 indicators: " + "; ".join(f"{cn}: {META[cn]['indicator']}" for cn in S2))
M_full = (B[train_pairs].astype(np.float64).T @ B[train_pairs].astype(np.float64)) / len(train_pairs)
if A.arm == "nn":
    torch.manual_seed(0)
    NET = torch.nn.Sequential(torch.nn.Linear(B.shape[1], 64), torch.nn.GELU(), torch.nn.Linear(64, 32), torch.nn.GELU(), torch.nn.Linear(32, 1)).to(DEV)
    torch.nn.init.zeros_(NET[-1].weight); torch.nn.init.zeros_(NET[-1].bias)
    NPAR = sum(p.numel() for p in NET.parameters()); log(f"nn params {NPAR}")
    Bt_gpu = torch.tensor(B, device=DEV)
    def nn_pack():   return torch.cat([p.detach().reshape(-1) for p in NET.parameters()]).cpu().numpy().astype(np.float64)
    def nn_unpack(v):
        i = 0
        for p in NET.parameters():
            n_ = p.numel(); p.data.copy_(torch.tensor(v[i:i+n_], dtype=p.dtype, device=DEV).reshape(p.shape)); i += n_
    def nn_f(): return NET(Bt_gpu).squeeze(1)
    Pfull = np.eye(NPAR); A0_NN = nn_pack()
NPARAM = NPAR if A.arm == "nn" else B.shape[1]
PW = (1.0 / np.maximum(col("n_genes"), 1.0)).astype(np.float32)
log(f"pair weight 1/n: mean {PW.mean():.3f}, n_genes>1 share {(PW < 1).mean()*100:.1f}%")

NIT = [0]; FROZEN = [None]; FIT_GI = [None]
def f_of(a, dt=np.float32):
    if A.arm == "nn":
        nn_unpack(a)
        with torch.no_grad(): return nn_f().cpu().numpy().astype(dt)
    return (B.astype(dt) @ a.astype(dt))
def phi_of(a, dt=np.float32): return np.exp(np.clip(f_of(a, dt), -CL, CL)).astype(dt)
def data_loss(a, want_grad=True):
    dt = np.float64 if DT[0] == torch.float64 else np.float32
    eta = f_of(a, dt); inside = (np.abs(eta) < CL); phi = np.exp(np.clip(eta, -CL, CL)).astype(dt)
    S = burden((phi * PW).astype(dt)).astype(dt); loss, gS = head(S, FIT_GI[0], FROZEN[0], want_grad=want_grad)
    if not want_grad: return loss, None
    gm = grad_phi(gS) * PW * phi * inside
    if A.arm == "nn":
        NET.zero_grad(); ft = nn_f(); ft.backward(gradient=torch.tensor(gm, device=DEV, dtype=ft.dtype))
        ga = torch.cat([p.grad.reshape(-1) for p in NET.parameters()]).cpu().numpy().astype(np.float64)
    else:
        ga = (B.astype(dt).T @ gm).astype(np.float64)
    return float(loss), ga
def objective(a, lam):
    loss, ga = data_loss(a)
    if not np.isfinite(loss): log(f"  non-finite loss |a|={np.abs(a).max():.2f} → 1e6"); return 1e6, (2 * lam * (Pfull @ a)).astype(np.float64)
    NIT[0] += 1
    pen = float(lam * (a @ (Pfull @ a)))
    if NIT[0] % 5 == 1: log(f"  it {NIT[0]} L {loss:.6f} pen {pen:.2e} |a|max={np.abs(a).max():.3f}")
    return loss + pen, ga + 2 * lam * (Pfull @ a)
def cfg_hash(gi, lam):
    h = hashlib.sha256(); h.update(CODE_SHA.encode()); h.update(json.dumps(dict(arm=A.arm, fold=A.fold, smoke=A.smoke, lam=lam, cl=CL, nknot=NKNOT, perm=A.perm_r,
                                                                             cols=COEF_NAMES, genes=[int(g) for g in gi]), sort_keys=True).encode())
    return h.hexdigest()[:16]
def fit(lam, gi, maxiter, tag):
    FIT_GI[0] = gi; FROZEN[0] = freeze_head(burden(PW), gi); NIT[0] = 0
    a = (A0_NN.copy() if A.arm == "nn" else np.zeros(NPARAM)).astype(np.float64)
    ck = f"{A.out}/{tag}.ckpt.npz"; H = cfg_hash(gi, lam)
    if os.path.exists(ck):
        z = np.load(ck, allow_pickle=True)
        if str(z["cfg"]) == H: a = z["a"].astype(np.float64); NIT[0] = int(z["nit"]); log(f"  RESUME {tag} nit {NIT[0]} |a|max={np.abs(a).max():.3f}")
        else: log(f"  ckpt ignored (cfg {str(z['cfg'])} != {H})")
    def cb(xk):
        if NIT[0] % 5 == 0: np.savez(ck + ".tmp.npz", a=xk, nit=NIT[0], cfg=H); os.replace(ck + ".tmp.npz", ck)
    res = minimize(objective, a, args=(lam,), jac=True, method="L-BFGS-B", callback=cb,
                   options=dict(maxiter=maxiter, maxfun=3 * maxiter, ftol=1e-9, gtol=1e-10, maxcor=20))
    a = res.x.astype(np.float64); log(f"  fit[{tag}] λ={lam:.4g} nit {res.nit} nfev {res.nfev} L+pen {res.fun:.6f} ok={res.success} {res.message}")
    os.path.exists(ck) and os.remove(ck)
    return a, res
def valid_J(a, gi_fit, gi_va):
    fr = freeze_head(burden(PW), gi_fit)
    S = burden(phi_of(a) * PW); J, _ = head(S, gi_va, fr, want_grad=False); return float(J)

t1 = time.time(); S_flat = burden(PW); log(f"burden(flat·1/n) {time.time()-t1:.1f}s")
FIT_GI[0] = train_gi; FROZEN[0] = freeze_head(S_flat, train_gi); log(f"frozen (π,τ²) train flat: { {t: (round(v[0],4), round(v[1],3)) for t, v in FROZEN[0].items()} }")
if A.smoke:
    rng0 = np.random.default_rng(7); a_t = (A0_NN.copy() if A.arm == "nn" else np.zeros(NPARAM)) + rng0.normal(0, 0.05, NPARAM); d = rng0.normal(0, 1, NPARAM); d /= np.linalg.norm(d)
    DT[0] = torch.float64
    if A.arm == "nn": NET.double(); Bt_gpu = Bt_gpu.double()
    lam_t = 1.0; f0, g0 = objective(a_t, lam_t); eps = 1e-3; f1, _ = objective(a_t + eps * d, lam_t); f2, _ = objective(a_t - eps * d, lam_t)
    DT[0] = torch.float32
    if A.arm == "nn": NET.float(); Bt_gpu = Bt_gpu.float()
    fd, an = (f1 - f2) / (2 * eps), float(g0 @ d); log(f"GRADCHECK fd={fd:.6e} analytic={an:.6e} rel={abs(fd-an)/(abs(fd)+1e-12):.3e}")
    require(abs(fd - an) / (abs(fd) + 1e-12) < 0.05, "gradient check failed")
    if A.arm == "spline":
        for b in BLOCKS:
            if b[1] == "spl": ev = np.linalg.eigvalsh(b[3]); require(ev.min() > 1e-10, f"P not PD for {b[0]}: min ev {ev.min():.2e}")
        log("P blocks PD ok")
    NIT[0] = 0
    if A.gradcheck_only: print("GRADCHECK_ONLY_DONE"); sys.exit(0)

a0 = (A0_NN.copy() if A.arm == "nn" else np.zeros(NPARAM)).astype(np.float64)
FIT_GI[0] = inner_tr; FROZEN[0] = freeze_head(S_flat, inner_tr)
L0, g0 = data_loss(a0)
if A.arm == "nn":
    idx = torch.tensor(train_pairs[::max(1, len(train_pairs) // 200000)], device=DEV)
    with torch.no_grad(): H = NET[:-1](Bt_gpu[idx]); H = torch.cat([H, torch.ones((H.shape[0], 1), device=DEV, dtype=H.dtype)], 1)
    Mh = (H.T @ H / H.shape[0]).cpu().numpy().astype(np.float64); gl_ = g0[-33:]; del H
    lam_strong = float(np.sqrt(gl_ @ Mh @ gl_) / (2 * EPS_F))
else:
    Pinv_g = np.linalg.solve(Pfull, g0); lam_strong = float(np.sqrt(Pinv_g @ M_full @ Pinv_g) / (2 * EPS_F))
log(f"λ강 = {lam_strong:.4g}  (|g0|={np.linalg.norm(g0):.3e}, ε={EPS_F}, L0={L0:.6f})")

TAG = A.arm + "_v7_" + ("smoke" if A.smoke else f"fold{A.fold}") + (f"_chr{A.smoke}" if A.smoke else "") + (f"_perm{A.perm_r}" if A.perm_r else "")
inner = {}
if A.lam is None:
    J_flat = valid_J(a0, inner_tr, inner_va); inner["inf"] = dict(J=J_flat)
    grid = [lam_strong / (10 ** k) for k in range(A.lam_grid)]
    k = 0
    while k < len(grid):
        lam = grid[k]; a_k, res_k = fit(lam, inner_tr, A.inner_maxiter, f"{TAG}_inner{k}")
        Jk = valid_J(a_k, inner_tr, inner_va); inner[f"{lam:.4g}"] = dict(J=Jk, nit=int(res_k.nit), amax=float(np.abs(a_k).max()), rms_f=float(np.sqrt((f_of(a_k)[train_pairs] ** 2).mean())))
        log(f"  inner λ={lam:.4g}: J_va {Jk:.6f} (flat {J_flat:.6f}) rms_f {inner[f'{lam:.4g}']['rms_f']:.4f}")
        if k == len(grid) - 1 and k == A.lam_grid - 1 and Jk < min(v["J"] for kk, v in inner.items() if kk != f"{lam:.4g}"):
            grid.append(lam / 10); log("  edge rule: weakest wins → extend grid ×1/10 once")
        k += 1
    cands = [(float("inf"), J_flat)] + [(lam, inner[f"{lam:.4g}"]["J"]) for lam in grid]
    cands.sort(key=lambda z: -z[0] if z[0] != float("inf") else -1e300)
    best = min(cands, key=lambda z: z[1])
    for lam, J in cands:
        if J - best[1] < TIE * abs(J_flat): chosen = lam; break
    log(f"  λ chosen = {chosen}  (best J {best[1]:.6f} at λ={best[0]})")
else:
    chosen = A.lam
if chosen == float("inf"):
    a = a0; res = None; log("  chosen ∞ → learned = flat")
else:
    a, res = fit(chosen, train_gi, A.maxiter, TAG)
phi = phi_of(a); S_learn = burden(phi * PW)
ho_l, ho_f = holdout(S_learn, test_gi), holdout(S_flat, test_gi)
ratio = {t: ho_l[t] / ho_f[t] for t in TRAITS}
infl = {t: float(np.median(gene_T(S_learn, test_gi, t) ** 2) / np.median(gene_T(S_flat, test_gi, t) ** 2)) for t in TRAITS}
eta = f_of(a); clamp_frac = float((np.abs(eta) >= CL).mean())
corr_diag = {cn: float(np.corrcoef(np.nan_to_num(col(cn)), np.log(phi))[0, 1]) if phi.std() > 0 else 0.0 for cn in DIAG}
null = {}
if A.null_perm:
    rng = np.random.default_rng(123); Zl = {t: zrow(S_learn[np.ix_(test_gi, TR[t]["cols"])]) for t in TRAITS}; Zf = {t: zrow(S_flat[np.ix_(test_gi, TR[t]["cols"])]) for t in TRAITS}
    rat = []
    for b in range(A.null_perm):
        perm = rng.permutation(NS); rr = []
        for t in TRAITS:
            full = FULL[t][perm]; m = ~np.isnan(FULL[t]); r_ = full[m]
            r_ = np.nan_to_num(r_ - np.nanmean(r_)); rr.append(np.abs(Zl[t] @ r_).mean() / (np.abs(Zf[t] @ r_).mean() + 1e-12))
        rat.append(rr)
    rat = np.array(rat); obs = np.array([ratio[t] for t in TRAITS]); pm = rat.mean(1)
    null = dict(B=A.null_perm, q95={t: float(np.quantile(rat[:, i], .95)) for i, t in enumerate(TRAITS)}, q95_mean=float(np.quantile(pm, .95)),
                p={t: float(((rat[:, i] >= obs[i]).sum() + 1) / (A.null_perm + 1)) for i, t in enumerate(TRAITS)}, p_mean=float(((pm >= obs.mean()).sum() + 1) / (A.null_perm + 1)))
out = dict(version="v7", code_sha256=CODE_SHA, arm=A.arm, mode="smoke" if A.smoke else f"fold{A.fold}", chr=sorted(use_chr, key=int), pairs=int(len(keys)), genes=NG,
           train_genes=int(len(train_gi)), inner_tr=int(len(inner_tr)), inner_va=int(len(inner_va)), test_genes=int(len(test_gi)), n_param=int(NPARAM), blocks=META,
           pair_weight="1/n_genes", penalty=dict(form="lambda a'Pa; spl P=G+Om/q, bin/ind P=1, lin P=E[x^2], nn P=I", CL=CL, EPS_F=EPS_F, TIE=TIE, lam_strong=lam_strong, lam_chosen=(None if chosen == float("inf") else chosen), inner=inner),
           frozen_head={t: [float(v[0]), float(v[1])] for t, v in FROZEN[0].items()}, perm_r=A.perm_r,
           nit=(int(res.nit) if res is not None else 0), converged=(bool(res.success) if res is not None else True), loss=(float(res.fun) if res is not None else None),
           holdout_learned=ho_l, holdout_flat=ho_f, ratio=ratio, ratio_mean=float(np.mean(list(ratio.values()))), inflation_medT2=infl,
           phi_q=[float(q) for q in np.quantile(phi, [0, .01, .5, .99, 1])], clamp_frac=clamp_frac, corr_logphi_diag=corr_diag, null_fixed_phi=null, sec=round(time.time() - T0))
np.save(f"{A.out}/{TAG}.a.npy", a); np.save(f"{A.out}/{TAG}.phi.npy", phi); json.dump(out, open(f"{A.out}/{TAG}.json", "w"), indent=1)
json.dump(dict(coef_names=COEF_NAMES), open(f"{A.out}/{TAG}.coef_names.json", "w"))
log("RESULT", json.dumps(out)); print("L1_DONE")
