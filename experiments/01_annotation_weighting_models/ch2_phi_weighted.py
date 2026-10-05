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


import csv, json, random, statistics as st, sys
random.seed(20260910)
D = _config_path("${PROJECT_ROOT}/work/run_l3b/out")
phi = {}
for line in open(f"{D}/G_phi.txt"):
    t = line.split()
    if len(t) >= 3 and t[1] == "weight":
        phi[t[0]] = [float(x) for x in t[2:]]
mphi = {g: st.mean(v) for g, v in phi.items()}
P = {}
for r in csv.DictReader(open(f"{D}/none.merged", encoding="utf-8", errors="replace"), delimiter="\t"):
    if r.get("Group") == "all":
        try: P[r["Region"]] = float(r["Pvalue"])
        except (ValueError, KeyError): pass
genes = sorted(set(P) & set(mphi))
G = len(genes); THR = 2.5e-6
assert G > 1400, f"조인 부족 G={G} (P {len(P)}, phi {len(mphi)})"
S = sum(mphi[g] for g in genes)
def count_pass(wmap):
    return sum(1 for g in genes if P[g] / wmap[g] < THR)
w_obs = {g: G * mphi[g] / S for g in genes}
n_unw = sum(1 for g in genes if P[g] < THR)
n_obs = count_pass(w_obs)
vals = [mphi[g] for g in genes]; null = []
for _ in range(200):
    random.shuffle(vals)
    wm = {g: G * v / S for g, v in zip(genes, vals)}
    null.append(count_pass(wm))
null_sorted = sorted(null)
p95 = null_sorted[int(0.95 * 200) - 1]
p_perm = (sum(1 for x in null if x >= n_obs) + 1) / 201
out = {"G": G, "threshold": THR, "n_unweighted": n_unw, "n_weighted": n_obs,
       "null_mean": st.mean(null), "null_sd": st.pstdev(null), "null_p95": p95, "null_max": max(null), "p_perm_ge": p_perm,
       "mean_phi_stats": {"min": min(vals), "median": st.median(vals), "max": max(vals), "sd": st.pstdev(vals)},
       "w_stats": {"min": min(w_obs.values()), "median": st.median(w_obs.values()), "max": max(w_obs.values())},
       "verdict_rule": "n_weighted - n_unweighted > null_p95 - n_unweighted → 전달 실증", "seed": 20260910, "B": 200}
tracked = ["ENSG00000129353.15","ENSG00000213892.12","ENSG00000186567.14",
           "ENSG00000130202.10","ENSG00000130204.13","ENSG00000104856.15","ENSG00000069399.15",
           "ENSG00000079805.19","ENSG00000142453.13","ENSG00000127616.22","ENSG00000129354.12"]
rows = []
for g in tracked:
    if g in P and g in w_obs:
        rows.append({"gene": g, "p": P[g], "mean_phi": mphi[g], "w": w_obs[g], "p_weighted": P[g] / w_obs[g],
                     "pass_unw": P[g] < THR, "pass_w": P[g] / w_obs[g] < THR})
out["tracked"] = rows
out["gained"] = sorted([g for g in genes if P[g] >= THR and P[g] / w_obs[g] < THR], key=lambda g: P[g])
out["lost"] = sorted([g for g in genes if P[g] < THR and P[g] / w_obs[g] >= THR], key=lambda g: P[g])
json.dump(out, open(f"{D}/ch2_phi_weighted.json", "w"), indent=1)
print(json.dumps({k: out[k] for k in ("G","n_unweighted","n_weighted","null_mean","null_sd","null_p95","null_max","p_perm_ge","mean_phi_stats","w_stats")}, ensure_ascii=False))
print("gained:", out["gained"]); print("lost:", out["lost"])
for r in rows: print(f"  {r['gene']:22s} p={r['p']:.3e} mphi={r['mean_phi']:.4f} w={r['w']:.3f} p'={r['p_weighted']:.3e} {'PASS' if r['pass_w'] else '-'}{'*' if r['pass_w'] and not r['pass_unw'] else ''}")
print("CH2_DONE")
