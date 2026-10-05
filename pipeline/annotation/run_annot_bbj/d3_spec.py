#!/usr/bin/env python3
import sys, json, gzip, os, math
phi, pipdir, traits, out = sys.argv[1:5]
m = json.load(open(phi)); assert m.get("version") == "v10" and list(m["coefficients"]) == ["shared"]
tr = [l.split("\t")[0] for l in open(traits).read().splitlines()[1:]]
num = den = 0.0; ncs = 0; nrows = 0
for t in tr:
    groups = {}
    with gzip.open(f"{pipdir}/{t}.tsv.gz", "rt") as fh:
        hdr = next(fh).rstrip("\n").split("\t"); ip, ic, ir = hdr.index("pip"), hdr.index("cs_id"), hdr.index("region")
        for line in fh:
            a = line.rstrip("\n").split("\t")
            if a[ic] == "-1": continue
            groups.setdefault((a[ir], a[ic]), []).append(float(a[ip])); nrows += 1
    for k, v in groups.items():
        if len(v) < 3: continue
        mu = sum(v) / len(v); var = sum((x - mu) ** 2 for x in v) / (len(v) - 1)
        num += mu * (1 - mu); den += var; ncs += 1
assert den > 0 and ncs > 0
s = num / den - 1
assert s > 0, f"concentration fit nonpositive {s}"
spec = dict(phi_model=os.path.abspath(phi), intercepts={k: float(v) for k, v in m["trait_intercepts"].items()},
            label_law="beta_mean_sigmoid", concentration=float(s),
            provenance=dict(method="pooled method of moments over CS (cs_id!=-1, n>=3): s = sum m(1-m)/sum var - 1 (ddof=1)",
                            n_cs=ncs, n_cs_rows=nrows, sum_m1m=num, sum_var=den,
                            delta_definition="delta multiplies the fitted additive phi score basis@coefficients (delta=1 == exported fitted phi)"))
missing = [t for t in tr if t not in spec["intercepts"]]
assert not missing, f"traits without intercept: {len(missing)}"
json.dump(spec, open(out, "w"), indent=1)
print(json.dumps(dict(concentration=s, n_cs=ncs, n_rows=nrows, traits=len(tr))))
