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
import sys, glob, collections, time
import numpy as np
import torch, torch.nn as nn

R = _config_path("${PROJECT_ROOT}/work/ref")
DEV = "cuda" if torch.cuda.is_available() else "cpu"
FOLD0 = {"7","13","16","19","21"}
STICKERS = ["splice_acceptor_variant","splice_donor_variant","stop_gained","frameshift_variant",
            "stop_lost","start_lost","inframe_insertion","inframe_deletion","missense_variant",
            "protein_altering_variant","splice_region_variant"]
THR = 2.5e-6

sid, eta, y = [], [], []
with open(f"{R}/l0/tchl.offset.tsv") as fh:
    hdr = fh.readline()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        sid.append(f[0]); eta.append(float(f[1])); y.append(float(f[2]))
eta = np.array(eta, dtype=np.float32); y = np.array(y, dtype=np.float32)
sample_ix = {s: i for i, s in enumerate(sid)}
print(f"samples with TCHL: {len(sid):,}")

sig = {}
import re
part_re = re.compile(r"\.part\d{3}$")
for fp in glob.glob(f"{R}/saige_step2_v4/tchl.chr*.part*"):
    if not part_re.search(fp): continue
    ch = fp.split(".chr")[1].split(".")[0]
    with open(fp) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        gi, pi = h.index("Region"), h.index("Pvalue")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            try: p = float(f[pi])
            except (ValueError, IndexError): continue
            if p < THR: sig[f[gi].split(".")[0]] = ch
train_seeds = {g for g, ch in sig.items() if ch not in FOLD0}
test_seeds  = {g for g, ch in sig.items() if ch in FOLD0}
print(f"v4 sig genes: {len(sig)} | train-fold seeds: {len(train_seeds)} | test-fold: {len(test_seeds)}")
assert len(train_seeds) >= 2, "GATE FAIL: too few training seeds"

want = train_seeds | test_seeds
gv = collections.defaultdict(list)
with open(f"{R}/deductive/l0_variants.tsv") as fh:
    fh.readline()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if f[0] in want:
            gv[f[0]].append((f[1], np.array(f[2:], dtype=np.float32)))
print(f"seed genes with L0 coding variants: {len(gv)} / {len(want)}")
kept_train = [g for g in train_seeds if g in gv]
kept_test  = [g for g in test_seeds if g in gv]
print(f"  usable train seeds: {len(kept_train)}  test: {len(kept_test)}")
assert kept_train, "GATE FAIL: no train seed has coding variants"

need = {k for g in gv for k, _ in gv[g]}
byc = collections.defaultdict(set)
for k in need: byc[k.split(":")[0].replace("chr","")].add(k)
DS = {}
for ch, keys in byc.items():
    for line in open(f"{R}/l0/chr{ch}.ds.tsv"):
        p = line.rstrip("\n").split("\t")
        key = "chr" + p[0]
        if key in keys:
            DS[key] = np.array(p[1:], dtype=np.float32)
print(f"dosage rows loaded: {len(DS):,} / {len(need):,}")
assert len(DS) == len(need), f"GATE FAIL: missing dosage rows {len(need)-len(DS)}"

import subprocess, os
sf = f"{R}/l0/band_samples.txt"
if not os.path.exists(sf):
    subprocess.run(f"bcftools query -l {R}/band_vcf/chr21.band.vcf.gz > {sf}", shell=True)
vcf_samples = [l.strip() for l in open(sf)]
assert len(vcf_samples) == N_SAMPLES, f"GATE FAIL: vcf samples {len(vcf_samples)}"
col_ix = np.array([i for i, s in enumerate(vcf_samples) if s in sample_ix])
row_map = np.array([sample_ix[vcf_samples[i]] for i in col_ix])
print(f"aligned samples: {len(col_ix):,}")
assert len(col_ix) == len(sid), "GATE FAIL: sample alignment"

eta_a = eta[row_map]; y_a = y[row_map]

genes = kept_train + kept_test
gene_slices = []; X_rows = []; ds_rows = []
for g in genes:
    idx0 = len(X_rows)
    for k, sv in gv[g]:
        X_rows.append(sv); ds_rows.append(DS[k][col_ix])
    gene_slices.append((idx0, len(X_rows)))
X = torch.tensor(np.stack(X_rows), device=DEV)
DSm = torch.tensor(np.stack(ds_rows), device=DEV)
ETA = torch.tensor(eta_a, device=DEV); Y = torch.tensor(y_a, device=DEV)
n_train_g = len(kept_train)
print(f"tensor: variants {X.shape[0]:,} x samples {DSm.shape[1]:,} on {DEV}")

class Phi(nn.Module):
    def __init__(s):
        super().__init__()
        s.net = nn.Sequential(nn.Linear(11,20), nn.LeakyReLU(), nn.Linear(20,20), nn.LeakyReLU(), nn.Linear(20,1))
    def forward(s, x): return nn.functional.softplus(s.net(x)).squeeze(-1)

torch.manual_seed(0)
phi = Phi().to(DEV)
beta = nn.Parameter(torch.zeros(n_train_g, device=DEV))
b0 = nn.Parameter(torch.zeros(1, device=DEV))
opt = torch.optim.AdamW(list(phi.parameters()) + [beta, b0], lr=1e-3)
resid = Y - ETA

losses = []
t0 = time.time()
for ep in range(200):
    opt.zero_grad()
    m = phi(X)
    S = []
    for gi in range(n_train_g):
        a, b = gene_slices[gi]
        S.append((m[a:b, None] * DSm[a:b]).sum(0))
    S = torch.stack(S)
    S = (S - S.mean(1, keepdim=True)) / (S.std(1, keepdim=True) + 1e-8)
    pred = b0 + (beta[:, None] * S).sum(0)
    loss = ((resid - pred) ** 2).mean()
    loss.backward(); opt.step()
    losses.append(loss.item())
    if ep % 50 == 0: print(f"ep {ep}: loss {loss.item():.6f}", flush=True)
dt = time.time() - t0
print(f"train time: {dt:.1f}s ({dt/200*1000:.0f} ms/epoch)")

drop = losses[0] - losses[-1]
mono = losses[-1] < losses[0] and losses[-1] <= min(losses) * 1.001
m_final = phi(X).detach()
m_std = m_final.std().item()
print(f"loss {losses[0]:.6f} -> {losses[-1]:.6f} (drop {drop:.6f})")
print(f"m_v: mean {m_final.mean():.4f} std {m_std:.4f} min {m_final.min():.4f} max {m_final.max():.4f}")
with torch.no_grad():
    rep = []
    for gi in range(n_train_g, len(genes)):
        a, b = gene_slices[gi]
        S_l = (m_final[a:b, None] * DSm[a:b]).sum(0)
        S_f = DSm[a:b].sum(0)
        for tag, S_ in (("learned", S_l), ("flat", S_f)):
            S_ = (S_ - S_.mean()) / (S_.std() + 1e-8)
            r = torch.dot(S_, resid - resid.mean()) / (len(resid) * resid.std() + 1e-8)
            rep.append((genes[gi], tag, r.item()))
print("holdout gene | corr(learned) vs corr(flat):")
for g in kept_test[:8]:
    rl = next(r for gg, t, r in rep if gg == g and t == "learned")
    rf = next(r for gg, t, r in rep if gg == g and t == "flat")
    print(f"  {g}  learned {rl:+.4f}  flat {rf:+.4f}")
assert drop > 0 and m_std > 1e-4, "GATE FAIL: no learning signal"
print("L0_SMOKE_OK")
