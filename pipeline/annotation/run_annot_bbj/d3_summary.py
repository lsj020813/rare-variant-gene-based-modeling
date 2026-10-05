#!/usr/bin/env python3
import json, sys, os, math
D = sys.argv[1]
def binom_ci(k, n, alpha=0.05):
    from math import lgamma, exp, log
    def betainc_cdf(x, a, b, steps=20000):
        if x <= 0: return 0.0
        if x >= 1: return 1.0
        s = 0.0; h = x / steps; c = exp(lgamma(a + b) - lgamma(a) - lgamma(b))
        for i in range(steps):
            t = (i + 0.5) * h; s += t ** (a - 1) * (1 - t) ** (b - 1)
        return min(1.0, c * s * h)
    def q(p, a, b):
        lo, hi = 0.0, 1.0
        for _ in range(60):
            mid = (lo + hi) / 2
            if betainc_cdf(mid, a, b) < p: lo = mid
            else: hi = mid
        return (lo + hi) / 2
    lo = 0.0 if k == 0 else q(alpha / 2, k, n - k + 1)
    hi = 1.0 if k == n else q(1 - alpha / 2, k + 1, n - k)
    return lo, hi
sim = None; status = "complete"
if os.path.exists(f"{D}/RESULT.json"):
    r = json.load(open(f"{D}/RESULT.json")); sim = r.get("simulation")
if sim is None and os.path.exists(f"{D}/simulation.partial.json"):
    sim = json.load(open(f"{D}/simulation.partial.json")); status = "partial"
assert sim, "no simulation output"
alpha = 0.05; power_target = 0.8; type1_max = 0.10
rows = []
def add(mode, delta, reps, note=""):
    n = len(reps); k = sum(1 for x in reps if x.get("significant"))
    lo, hi = binom_ci(k, n, alpha) if n else (None, None)
    rec = [x.get("phi_recovery_correlation") for x in reps if x.get("phi_recovery_correlation") is not None]
    rows.append(dict(mode=mode, delta=delta, n_reps=n, rejections=k, rejection_rate=(k / n if n else None), ci95_lo=lo, ci95_hi=hi,
                     power_ok=(lo is not None and lo >= power_target) if mode == "power" else None,
                     type1_upper_ok=(hi is not None and hi <= type1_max) if mode == "null" else None,
                     mean_phi_recovery_corr=(sum(rec) / len(rec) if rec else None), note=note))
for mode in ("null", "oracle", "power"):
    for e in sim.get(mode, []):
        add(mode, e.get("delta"), e.get("rows", []))
cur = sim.get("current")
if cur and cur.get("rows"):
    add(cur["mode"], cur.get("delta"), cur["rows"], note="in progress")
elig = [r["delta"] for r in rows if r["mode"] == "power" and r["power_ok"] and r["note"] == ""]
meta = dict(status=status, delta_star_candidate=(min(elig) if elig else None), delta_definition="delta multiplies the fitted primary_f0 phi additive score (basis@coef; delta=1 == fitted phi; ||coef||_2 = 88.34)",
            null_upper95=[r["ci95_hi"] for r in rows if r["mode"] == "null"], type1_max=type1_max, power_target=power_target, alpha=alpha,
            conservative_note="training/evaluation restricted to the 4 primary test chromosomes (7,13,16,19): reduced training and test size -> power estimate is conservative")
keys = list(rows[0].keys()) if rows else []
with open(f"{D}/d3_summary.tsv", "w") as o:
    o.write("\t".join(keys) + "\n")
    for r in rows: o.write("\t".join("" if r[k] is None else (f"{r[k]:.4g}" if isinstance(r[k], float) else str(r[k])) for k in keys) + "\n")
json.dump(meta, open(f"{D}/d3_summary.meta.json", "w"), indent=1)
print(json.dumps(meta)); print("\n".join("\t".join(str(r[k]) for k in keys) for r in rows))
