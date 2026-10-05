#!/usr/bin/env python3
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import glob, re, collections, time
import numpy as np
import torch, torch.nn as nn

R = _config_path("${PROJECT_ROOT}/work/ref")
DEV = "cuda" if torch.cuda.is_available() else "cpu"
FOLD0 = {"7","13","16","19","21"}
STICKERS = ["splice_acceptor_variant","splice_donor_variant","stop_gained","frameshift_variant",
            "stop_lost","start_lost","inframe_insertion","inframe_deletion","missense_variant",
            "protein_altering_variant","splice_region_variant"]
TRAITS = ["tchl","htn","dm","lip"]

gv = collections.defaultdict(list)
with open(f"{R}/deductive/l0_variants.tsv") as fh:
    fh.readline()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        gv[f[0]].append((f[1], np.array(f[2:], dtype=np.float32)))
gchr = {g: vs[0][0].split(":")[0].replace("chr","") for g, vs in gv.items()}
train_genes = sorted(g for g in gv if gchr[g] not in FOLD0)
test_genes  = sorted(g for g in gv if gchr[g] in FOLD0)
print(f"genes: train {len(train_genes):,}  test(fold0) {len(test_genes):,}")

need = {k for g in gv for k, _ in gv[g]}
byc = collections.defaultdict(set)
for k in need: byc[k.split(":")[0].replace("chr","")].add(k)
DS = {}
for ch, keys in byc.items():
    for line in open(f"{R}/l0/chr{ch}.ds.tsv"):
        p = line.rstrip("\n").split("\t")
        key = "chr" + p[0]
        if key in keys and key not in DS:
            DS[key] = np.frombuffer(np.array(p[1:], dtype=np.float32).tobytes(), dtype=np.float32)
print(f"dosage rows: {len(DS):,}/{len(need):,}")
assert len(DS) == len(need), "GATE FAIL: dosage missing"

vcf_samples = [l.strip() for l in open(f"{R}/l0/band_samples.txt")]
tr = {}
for t in TRAITS:
    sid, eta, y = [], [], []
    with open(f"{R}/l0/{t}.offset.tsv") as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            sid.append(f[0]); eta.append(float(f[1])); y.append(float(f[2]))
    ix = {s: i for i, s in enumerate(sid)}
    cols = np.array([i for i, s in enumerate(vcf_samples) if s in ix])
    rows = np.array([ix[vcf_samples[i]] for i in cols])
    tr[t] = dict(cols=cols, resid=np.array(y, np.float32)[rows] - np.array(eta, np.float32)[rows])
    print(f"{t}: n={len(cols):,}")

genes = train_genes + test_genes
X_rows, ds_rows, slices = [], [], []
for g in genes:
    a = len(X_rows)
    for k, sv in gv[g]:
        X_rows.append(sv); ds_rows.append(DS[k])
    slices.append((a, len(X_rows)))
X = torch.tensor(np.stack(X_rows), device=DEV)
DSfull = torch.tensor(np.stack(ds_rows), device=DEV)
n_tr_g = len(train_genes)
print(f"tensor: {X.shape[0]:,} variants x configured samples on {DEV}")

def burden(m, cols):
    D = DSfull[:, cols]
    return torch.stack([ (m[a:b, None] * D[a:b]).sum(0) for a, b in slices ])

def fit_and_eval(phi_fn, params, tag, epochs=150, lam=10.0):
    opt = torch.optim.AdamW(params, lr=5e-3)
    losses = []
    t0 = time.time()
    for ep in range(epochs):
        opt.zero_grad()
        m = phi_fn()
        loss = 0.0
        for t in TRAITS:
            cols = torch.tensor(tr[t]["cols"], device=DEV)
            resid = torch.tensor(tr[t]["resid"], device=DEV)
            S = burden(m, cols)[:n_tr_g]
            S = (S - S.mean(1, keepdim=True)) / (S.std(1, keepdim=True) + 1e-8)
            G = S @ S.T + lam * torch.eye(n_tr_g, device=DEV)
            beta = torch.linalg.solve(G, S @ resid)
            pred = beta @ S
            loss = loss + ((resid - pred) ** 2).mean()
        loss.backward(); opt.step()
        losses.append(loss.item())
        if ep % 50 == 0: print(f"[{tag}] ep {ep}: {loss.item():.5f}", flush=True)
    print(f"[{tag}] {time.time()-t0:.0f}s  loss {losses[0]:.5f} -> {losses[-1]:.5f}")
    with torch.no_grad():
        m = phi_fn()
        cols = torch.tensor(tr["tchl"]["cols"], device=DEV)
        resid = torch.tensor(tr["tchl"]["resid"], device=DEV)
        rs_l, rs_f = [], []
        D = DSfull[:, cols]
        for gi in range(n_tr_g, len(genes)):
            a, b = slices[gi]
            for S_, acc in (((m[a:b, None]*D[a:b]).sum(0), rs_l), (D[a:b].sum(0), rs_f)):
                S_ = (S_ - S_.mean()) / (S_.std() + 1e-8)
                acc.append((torch.dot(S_, resid - resid.mean()) / (len(resid) * resid.std() + 1e-8)).abs().item())
        rl, rf = float(np.mean(rs_l)), float(np.mean(rs_f))
        print(f"[{tag}] holdout mean|corr| learned {rl:.5f} vs flat {rf:.5f}  ratio {rl/max(rf,1e-12):.3f}")
    return losses, m

alpha = nn.Parameter(torch.zeros(11, device=DEV))
losses_a, m_a = fit_and_eval(lambda: torch.exp(X @ alpha), [alpha], "loglinear")
print("alpha:", {s: round(a,3) for s, a in zip(STICKERS, alpha.detach().cpu().tolist())})

mlp = nn.Sequential(nn.Linear(11,20), nn.LeakyReLU(), nn.Linear(20,20), nn.LeakyReLU(), nn.Linear(20,1)).to(DEV)
losses_b, m_b = fit_and_eval(lambda: nn.functional.softplus(mlp(X)).squeeze(-1), list(mlp.parameters()), "mlp")

assert losses_a[-1] < losses_a[0], "GATE FAIL: loglinear no descent"
assert alpha.detach().abs().max().item() > 1e-3, "GATE FAIL: alpha constant"
print("L0_SMOKE_V2_OK")
