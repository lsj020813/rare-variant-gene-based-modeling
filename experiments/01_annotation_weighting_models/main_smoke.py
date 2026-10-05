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
import collections, time, json, itertools
import numpy as np
import torch, torch.nn as nn

R = _config_path("${PROJECT_ROOT}/work/ref")
DEV = "cuda"
FOLDS = {0:{"7","13","16","19","21"}, 1:{"2","8","9","10","15"}, 2:{"1","6","14"},
         3:{"3","12","17","22"}, 4:{"4","5","11","18","20"}}
STICKERS = ["splice_acceptor_variant","splice_donor_variant","stop_gained","frameshift_variant",
            "stop_lost","start_lost","inframe_insertion","inframe_deletion","missense_variant",
            "protein_altering_variant","splice_region_variant"]
TRAITS = ["tchl","htn","dm","lip"]
EPOCHS, LAMS = 3, (10.0,)

CACHE = f"{R}/l0/l0_tensors.pt"
gv = collections.defaultdict(list)
with open(f"{R}/deductive/l0_variants.tsv") as fh:
    fh.readline()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        gv[f[0]].append((f[1], np.array(f[2:], dtype=np.float32)))
gchr = {g: vs[0][0].split(":")[0].replace("chr","") for g, vs in gv.items()}
genes_all = sorted(gv)

need = {k for g in gv for k, _ in gv[g]}
byc = collections.defaultdict(set)
for k in need: byc[k.split(":")[0].replace("chr","")].add(k)
import os
USE_CACHE = os.path.exists(CACHE)
DS = {}
for ch, keys in (byc.items() if not USE_CACHE else []):
    for line in open(f"{R}/l0/chr{ch}.ds.tsv"):
        p = line.rstrip("\n").split("\t")
        key = "chr" + p[0]
        if key in keys and key not in DS:
            DS[key] = np.array(p[1:], dtype=np.float16)
if not USE_CACHE:
    assert len(DS) == len(need), f"GATE FAIL dosage {len(DS)}/{len(need)}"
    print(f"dosage: {len(DS):,} variants x configured samples fp16", flush=True)
else:
    print("loading tensor cache", flush=True)

vcf_samples = [l.strip() for l in open(f"{R}/l0/band_samples.txt")]
tr = {}
for t in TRAITS:
    ix, eta, y = {}, [], []
    with open(f"{R}/l0/{t}.offset.tsv") as fh:
        fh.readline()
        for n, line in enumerate(fh):
            f = line.rstrip("\n").split("\t")
            ix[f[0]] = n; eta.append(float(f[1])); y.append(float(f[2]))
    eta = np.array(eta, np.float32); y = np.array(y, np.float32)
    cols = np.array([j for j, s in enumerate(vcf_samples) if s in ix])
    rows = np.array([ix[vcf_samples[j]] for j in cols])
    tr[t] = dict(cols=torch.tensor(cols, device=DEV),
                 resid=torch.tensor(y[rows]-eta[rows], device=DEV))
    print(f"{t}: n={len(cols):,}", flush=True)

if USE_CACHE:
    _c = torch.load(CACHE)
    X, Dh, GID = _c["X"].to(DEV), _c["Dh"].to(DEV), _c["GID"].to(DEV)
X_rows, ds_rows, gene_id = [], [], []
for gi, g in enumerate(genes_all):
    if USE_CACHE: break
    for k, sv in gv[g]:
        X_rows.append(sv); ds_rows.append(DS[k]); gene_id.append(gi)
if not USE_CACHE:
    X = torch.tensor(np.stack(X_rows), device=DEV)
    Dh = torch.tensor(np.stack(ds_rows), device=DEV)
    GID = torch.tensor(gene_id, device=DEV)
    torch.save({"X": X.cpu(), "Dh": Dh.cpu(), "GID": GID.cpu()}, CACHE)
    print("tensor cache saved", flush=True)
NG, NS = len(genes_all), N_SAMPLES
print(f"tensor: {X.shape[0]:,} x {NS:,} on {DEV} (fp16 D)", flush=True)

BLK = 2048
def burden_nograd(m, cols=None):
    with torch.no_grad():
        ns = NS if cols is None else len(cols)
        S = torch.zeros(NG, ns, device=DEV)
        for a in range(0, Dh.shape[0], BLK):
            b = min(a + BLK, Dh.shape[0])
            Db = Dh[a:b] if cols is None else Dh[a:b][:, cols]
            S.index_add_(0, GID[a:b], m[a:b, None] * Db.float())
        return S

def grad_m_from_gS(gS, cols=None):
    with torch.no_grad():
        gm = torch.zeros(Dh.shape[0], device=DEV)
        for a in range(0, Dh.shape[0], BLK):
            b = min(a + BLK, Dh.shape[0])
            Db = Dh[a:b] if cols is None else Dh[a:b][:, cols]
            gm[a:b] = (gS[GID[a:b]] * Db.float()).sum(1)
        return gm

def zrow(S):
    return (S - S.mean(1, keepdim=True)) / (S.std(1, keepdim=True) + 1e-8)

def holdout_score(m, test_gi):
    with torch.no_grad():
        S_all = burden_nograd(m)
        out = {}
        for t in TRAITS:
            cols, resid = tr[t]["cols"], tr[t]["resid"]
            rn = resid - resid.mean()
            S = zrow(S_all[test_gi][:, cols])
            out[t] = (S @ rn).abs().mean().item() / (len(resid) * resid.std().item() + 1e-8)
        return out

def train_arm(phi_fn, params, train_gi, lam, tag, fold):
    opt = torch.optim.AdamW(params, lr=5e-3)
    for ep in range(EPOCHS):
        opt.zero_grad()
        for t in TRAITS:
            cols_full, resid_full = tr[t]["cols"], tr[t]["resid"]
            bidx = torch.randperm(len(cols_full), device=DEV)[:30000]
            cols, resid = cols_full[bidx], resid_full[bidx]
            m = phi_fn()
            S0 = burden_nograd(m.detach(), cols)[train_gi]
            S_leaf = S0.detach().requires_grad_(True)
            S = zrow(S_leaf)
            G = S @ S.T + lam * torch.eye(len(train_gi), device=DEV)
            beta = torch.linalg.solve(G, S @ resid)
            loss = ((resid - beta @ S) ** 2).mean()
            loss.backward()
            gS_full = torch.zeros(NG, len(cols), device=DEV)
            gS_full[train_gi] = S_leaf.grad
            gm = grad_m_from_gS(gS_full, cols)
            m.backward(gradient=gm)
            del m, S0, S_leaf, S, G, beta, loss, gS_full, gm
            torch.cuda.empty_cache()
        opt.step()
        if ep % 100 == 0: print(f"[f{fold} {tag} lam{lam}] ep {ep}", flush=True)

res = {}
for fold in range(1):
    test_chr = FOLDS[fold]
    test_gi  = torch.tensor([i for i,g in enumerate(genes_all) if gchr[g] in test_chr], device=DEV)
    tr_genes = [i for i,g in enumerate(genes_all) if gchr[g] not in test_chr]
    rng = np.random.default_rng(fold)
    perm = rng.permutation(len(tr_genes))
    inner_tr = torch.tensor([tr_genes[i] for i in perm[:int(0.8*len(tr_genes))]], device=DEV)
    inner_va = torch.tensor([tr_genes[i] for i in perm[int(0.8*len(tr_genes)):]], device=DEV)
    full_tr  = torch.tensor(tr_genes, device=DEV)
    print(f"=== fold {fold}: train {len(tr_genes):,} test {len(test_gi):,} ===", flush=True)

    for tag, make in (("loglinear", lambda: ([nn.Parameter(torch.zeros(11, device=DEV))],)),
                      ("mlp",       lambda: (None,))):
        if tag == "loglinear":
            best_lam, best_s = None, -1
            for lam in LAMS:
                a = nn.Parameter(torch.zeros(11, device=DEV))
                train_arm(lambda: torch.exp(X @ a), [a], inner_tr, lam, f"{tag}-sel", fold)
                s = np.mean(list(holdout_score(torch.exp(X @ a).detach(), inner_va).values()))
                print(f"[f{fold}] lam {lam}: inner {s:.6f}", flush=True)
                if s > best_s: best_s, best_lam = s, lam
            LAM = best_lam
            print(f"[f{fold}] selected lam {LAM}", flush=True)
            alpha = nn.Parameter(torch.zeros(11, device=DEV))
            train_arm(lambda: torch.exp(X @ alpha), [alpha], full_tr, LAM, tag, fold)
            m_fin = torch.exp(X @ alpha).detach()
            res[f"f{fold}.alpha"] = {s: round(v,4) for s,v in zip(STICKERS, alpha.detach().cpu().tolist())}
        else:
            mlp = nn.Sequential(nn.Linear(11,20), nn.LeakyReLU(), nn.Linear(20,20),
                                nn.LeakyReLU(), nn.Linear(20,1)).to(DEV)
            train_arm(lambda: nn.functional.softplus(mlp(X)).squeeze(-1),
                      list(mlp.parameters()), full_tr, LAM, tag, fold)
            m_fin = nn.functional.softplus(mlp(X)).squeeze(-1).detach()
        sc_l = holdout_score(m_fin, test_gi)
        sc_f = holdout_score(torch.ones(X.shape[0], device=DEV), test_gi)
        res[f"f{fold}.{tag}"] = {t: {"learned": sc_l[t], "flat": sc_f[t]} for t in TRAITS}
        print(f"[f{fold} {tag}] " + " ".join(f"{t}:{sc_l[t]/max(sc_f[t],1e-12):.3f}" for t in TRAITS), flush=True)

json.dump(res, open(f"{R}/l0/l0_main_result.json","w"), indent=1)
print("L0_MAIN_COMPLETE")
