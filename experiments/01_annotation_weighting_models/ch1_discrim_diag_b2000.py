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


import numpy as np, csv, glob, json, math, sys
rng = np.random.default_rng(20260910)
W = _config_path("${PROJECT_ROOT}/work")
z = np.load(f"{W}/ref/annot/cache/fm_all.npz", allow_pickle=True)
m = z["chr"] == 19
X = z["X"][m].astype(np.float64); key = z["key37"][m]; gene = z["gene"][m]; cols = [str(c) for c in z["cols"]]
print("cache chr19:", X.shape, "genes", len(set(gene)))
sz = {}
for f in glob.glob(f"{W}/run_l3b/out/arm_none/part*.singleAssoc.txt"):
    for r in csv.DictReader(open(f), delimiter="\t"):
        try:
            v = float(r["var"]); 
            if v <= 0: continue
            sz[f'{r["CHR"]}:{r["POS"]}:{r["Allele1"]}:{r["Allele2"]}'] = abs(float(r["Tstat"]) / math.sqrt(v))
        except (ValueError, KeyError): pass
print("singleAssoc variants:", len(sz))
extra = {}
def load_tsv(path, names):
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        idx = [hdr.index(n) for n in names]
        for line in fh:
            t = line.rstrip("\n").split("\t")
            extra.setdefault(t[0], {}).update({n: t[i] for n, i in zip(names, idx)})
load_tsv(f"{W}/ref/annot/extract/chr19.div.tsv", ["div1","div2","div3","mutdens","cons_v2"])
try:
    with open(f"{W}/ref/annot/extract/chr19.cadd.tsv") as fh:
        hdr = fh.readline().rstrip("\n").split("\t"); ci = [i for i,h in enumerate(hdr) if "cadd" in h.lower()]
        ci = ci[0] if ci else 1
        for line in fh:
            t = line.rstrip("\n").split("\t"); extra.setdefault(t[0], {})["cadd"] = t[ci]
    have_cadd = True
except FileNotFoundError:
    have_cadd = False
print("extra keyed:", len(extra), "cadd:", have_cadd)
keep = np.array([k in sz for k in key]); print("matched |z|:", int(keep.sum()), "/", len(key))
X = X[keep]; key = key[keep]; gene = gene[keep]; y = np.array([sz[k] for k in key])
exclude = {"maf","r2","avg_cs","is_typed","is_indel","n_genes","gh_n_genes"}
names = [c for c in cols if c not in exclude]
M = {c: X[:, cols.index(c)] for c in names}
for n in ["div1","div2","div3","mutdens","cons_v2"] + (["cadd"] if have_cadd else []):
    M[n] = np.array([float(extra.get(k, {}).get(n) or "nan") if extra.get(k, {}).get(n, "") != "" else np.nan for k in key])
    names.append(n)
K = len(names); print("K columns tested:", K)
gidx = {}
for i, g in enumerate(gene): gidx.setdefault(g, []).append(i)
groups = [np.array(v) for v in gidx.values() if len(v) >= 3]
def within_rank(v):
    out = np.full(len(v), np.nan)
    for ix in groups:
        s = v[ix]; ok = ~np.isnan(s)
        if ok.sum() >= 3:
            r = np.empty(ok.sum()); r[np.argsort(s[ok], kind="stable")] = np.arange(ok.sum())
            out[ix[ok]] = r - r.mean()
    return out
ry = within_rank(y)
def corr(a, b):
    ok = ~np.isnan(a) & ~np.isnan(b); 
    if ok.sum() < 100: return np.nan, int(ok.sum())
    a, b = a[ok], b[ok]; a = a - a.mean(); b = b - b.mean()
    return float((a*b).sum() / math.sqrt((a*a).sum() * (b*b).sum())), int(ok.sum())
B = 2000; ynull = np.empty((B, len(y)))
for b in range(B):
    yy = y.copy()
    for ix in groups: yy[ix] = yy[rng.permutation(ix)]
    ynull[b] = within_rank(yy)
res = []
for n in names:
    rx = within_rank(M[n]); rho, nn = corr(rx, ry)
    if np.isnan(rho): res.append((n, np.nan, nn, np.nan, np.nan)); continue
    nul = np.array([corr(rx, ynull[b])[0] for b in range(B)])
    p = (np.sum(np.abs(nul) >= abs(rho)) + 1) / (B + 1)
    res.append((n, rho, nn, float(np.nanstd(nul)), p))
alpha = 0.05 / K
res.sort(key=lambda t: -abs(t[1]) if not np.isnan(t[1]) else 0)
print(f"  {'열':22s} {'rho(유전자내)':>13s} {'n':>7s} {'null_sd':>8s} {'p_perm':>7s}  판정(|rho|>=0.02 & p<{alpha:.4f})")
passed = []
for n, rho, nn, sd, p in res:
    ok = (not np.isnan(rho)) and abs(rho) >= 0.02 and p < alpha
    if ok: passed.append(n)
    print(f"  {n:22s} {rho:+13.4f} {nn:7d} {sd:8.4f} {p:7.4f}  {'★ 통과' if ok else ''}")
json.dump({"K": K, "alpha": alpha, "B": B, "n_variants": int(len(y)), "n_genes": len(groups),
           "rows": [{"col": n, "rho": None if np.isnan(r) else r, "n": nn, "null_sd": None if np.isnan(s) else s, "p_perm": None if np.isnan(p) else p} for n, r, nn, s, p in res],
           "passed": passed}, open(f"{W}/run_l3b/out/ch1_discrim_diag_B2000.json", "w"), indent=1)
print("PASSED:", passed); print("DIAG_DONE")
