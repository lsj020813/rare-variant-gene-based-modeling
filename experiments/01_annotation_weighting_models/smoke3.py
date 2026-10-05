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
import collections, time
import numpy as np
import torch, torch.nn as nn

R = _config_path("${PROJECT_ROOT}/work/ref")
DEV = "cuda" if torch.cuda.is_available() else "cpu"
FOLD0 = {"7","13","16","19","21"}
STICKERS = ["splice_acceptor_variant","splice_donor_variant","stop_gained","frameshift_variant",
            "stop_lost","start_lost","inframe_insertion","inframe_deletion","missense_variant",
            "protein_altering_variant","splice_region_variant"]
TRAITS = ["tchl","htn","dm","lip"]
NSUB = 20000

gv = collections.defaultdict(list)
with open(f"{R}/deductive/l0_variants.tsv") as fh:
    fh.readline()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        gv[f[0]].append((f[1], np.array(f[2:], dtype=np.float32)))
gchr = {g: vs[0][0].split(":")[0].replace("chr","") for g, vs in gv.items()}
train_genes = sorted(g for g in gv if gchr[g] not in FOLD0)
test_genes  = sorted(g for g in gv if gchr[g] in FOLD0)
genes = train_genes + test_genes
n_tr_g = len(train_genes)
print(f"genes: train {n_tr_g:,} test {len(test_genes):,}", flush=True)

need = {k for g in gv for k, _ in gv[g]}
byc = collections.defaultdict(set)
for k in need: byc[k.split(":")[0].replace("chr","")].add(k)

vcf_samples = [l.strip() for l in open(f"{R}/l0/band_samples.txt")]
rng = np.random.default_rng(0)
sub = np.sort(rng.choice(N_SAMPLES, NSUB, replace=False))
sub_set = set(sub.tolist())

DS = {}
for ch, keys in byc.items():
    for line in open(f"{R}/l0/chr{ch}.ds.tsv"):
        p = line.rstrip("\n").split("\t")
        key = "chr" + p[0]
        if key in keys and key not in DS:
            DS[key] = np.array(p[1:], dtype=np.float32)[sub]
print(f"dosage rows: {len(DS):,}/{len(need):,} (subsampled to {NSUB:,})", flush=True)
assert len(DS) == len(need)

tr = {}
sub_names = [vcf_samples[i] for i in sub]
for t in TRAITS:
    ix = {}
    eta, y = [], []
    with open(f"{R}/l0/{t}.offset.tsv") as fh:
        fh.readline()
        for n, line in enumerate(fh):
            f = line.rstrip("\n").split("\t")
            ix[f[0]] = n; eta.append(float(f[1])); y.append(float(f[2]))
    eta = np.array(eta, np.float32); y = np.array(y, np.float32)
    cols, rows = [], []
    for j, s in enumerate(sub_names):
        if s in ix:
            cols.append(j); rows.append(ix[s])
    cols = np.array(cols); rows = np.array(rows)
    tr[t] = dict(cols=torch.tensor(cols, device=DEV),
                 resid=torch.tensor(y[rows] - eta[rows], device=DEV))
    print(f"{t}: n={len(cols):,}", flush=True)

X_rows, ds_rows, gene_id = [], [], []
for gi, g in enumerate(genes):
    for k, sv in gv[g]:
        X_rows.append(sv); ds_rows.append(DS[k]); gene_id.append(gi)
X = torch.tensor(np.stack(X_rows), device=DEV)
D = torch.tensor(np.stack(ds_rows), device=DEV)
GID = torch.tensor(gene_id, device=DEV)
NG = len(genes)
print(f"tensor: {X.shape[0]:,} x {NSUB:,} on {DEV}", flush=True)

def burden(m):
    S = torch.zeros(NG, NSUB, device=DEV)
    S.index_add_(0, GID, m[:, None] * D)
    return S

def fit_and_eval(phi_fn, params, tag, epochs=150, lam=10.0):
    opt = torch.optim.AdamW(params, lr=5e-3)
    losses = []
    t0 = time.time()
    for ep in range(epochs):
        opt.zero_grad()
        m = phi_fn()
        S_all = burden(m)
        loss = 0.0
        for t in TRAITS:
            S = S_all[:n_tr_g][:, tr[t]["cols"]]
            S = (S - S.mean(1, keepdim=True)) / (S.std(1, keepdim=True) + 1e-8)
            resid = tr[t]["resid"]
            G = S @ S.T + lam * torch.eye(n_tr_g, device=DEV)
            beta = torch.linalg.solve(G, S @ resid)
            loss = loss + ((resid - beta @ S) ** 2).mean()
        loss.backward(); opt.step()
        losses.append(loss.item())
        if ep % 30 == 0: print(f"[{tag}] ep {ep}: {loss.item():.5f} ({time.time()-t0:.0f}s)", flush=True)
    print(f"[{tag}] total {time.time()-t0:.0f}s  loss {losses[0]:.5f} -> {losses[-1]:.5f}", flush=True)
    with torch.no_grad():
        m = phi_fn()
        S_all = burden(m)
        cols, resid = tr["tchl"]["cols"], tr["tchl"]["resid"]
        rn = resid - resid.mean()
        out = {}
        for name, Sx in (("learned", S_all), ("flat", burden(torch.ones_like(m)))):
            S = Sx[n_tr_g:][:, cols]
            S = (S - S.mean(1, keepdim=True)) / (S.std(1, keepdim=True) + 1e-8)
            r = (S @ rn) / (len(resid) * resid.std() + 1e-8)
            out[name] = r.abs().mean().item()
        print(f"[{tag}] holdout mean|corr|: learned {out['learned']:.5f} flat {out['flat']:.5f} ratio {out['learned']/max(out['flat'],1e-12):.3f}", flush=True)
    return losses

alpha = nn.Parameter(torch.zeros(11, device=DEV))
la = fit_and_eval(lambda: torch.exp(X @ alpha), [alpha], "loglinear")
print("alpha:", {s: round(a,3) for s, a in zip(STICKERS, alpha.detach().cpu().tolist())}, flush=True)

mlp = nn.Sequential(nn.Linear(11,20), nn.LeakyReLU(), nn.Linear(20,20), nn.LeakyReLU(), nn.Linear(20,1)).to(DEV)
lb = fit_and_eval(lambda: nn.functional.softplus(mlp(X)).squeeze(-1), list(mlp.parameters()), "mlp")

assert la[-1] < la[0], "GATE FAIL: no descent"
assert alpha.detach().abs().max().item() > 1e-3, "GATE FAIL: alpha constant"
print("L0_SMOKE_V3_OK")
