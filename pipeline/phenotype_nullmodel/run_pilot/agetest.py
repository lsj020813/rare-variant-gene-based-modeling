import os
import json

def required(name):
    value = os.environ.get(name, "")
    if not value.strip():
        raise RuntimeError("Required environment variable is missing or empty: " + name)
    return value

import csv, math
import statistics as st

def fit(X, y, iters=30):
    n, p = len(X), len(X[0])
    b = [0.0]*p
    for _ in range(iters):
        g = [0.0]*p
        H = [[0.0]*p for _ in range(p)]
        ll = 0.0
        for i in range(n):
            z = sum(b[k]*X[i][k] for k in range(p))
            z = max(-30, min(30, z))
            mu = 1/(1+math.exp(-z))
            w = max(mu*(1-mu), 1e-9)
            r = y[i]-mu
            ll += y[i]*math.log(max(mu,1e-12)) + (1-y[i])*math.log(max(1-mu,1e-12))
            for a in range(p):
                g[a] += X[i][a]*r
                for c2 in range(a, p): H[a][c2] += X[i][a]*X[i][c2]*w
        for a in range(p):
            for c2 in range(a): H[a][c2] = H[c2][a]
        M = [row[:]+[g[i]] for i, row in enumerate(H)]
        for c2 in range(p):
            piv = max(range(c2, p), key=lambda r2: abs(M[r2][c2]))
            M[c2], M[piv] = M[piv], M[c2]
            if abs(M[c2][c2]) < 1e-12: continue
            for r2 in range(p):
                if r2 == c2: continue
                f = M[r2][c2]/M[c2][c2]
                for k in range(c2, p+1): M[r2][k] -= f*M[c2][k]
        d = [M[i][p]/M[i][i] if abs(M[i][i])>1e-12 else 0.0 for i in range(p)]
        b = [b[i]+d[i] for i in range(p)]
        if max(abs(x) for x in d) < 1e-8: break
    return b, ll

for t in ('htn','dm','lip'):
    rows = list(csv.DictReader(open(os.path.join(required("PHENO_V2_DIR"), f"{t}_v2.tsv")), delimiter='\t'))
    y, A, S = [], [], []
    for r in rows:
        try:
            a = float(r['age']); s = float(r['sex_male']); yy = float(r['y'])
        except (ValueError, KeyError, TypeError):
            continue
        y.append(yy); A.append(a); S.append(s)
    m = sum(A)/len(A)
    ac = [(a-m) for a in A]
    X0 = [[1.0, ac[i], S[i]] for i in range(len(y))]
    X1 = [[1.0, ac[i], ac[i]**2, S[i]] for i in range(len(y))]
    b0, ll0 = fit(X0, y)
    b1, ll1 = fit(X1, y)
    lr = 2*(ll1-ll0)
    from math import erfc, sqrt
    p = erfc(sqrt(max(lr,0)/2))
    print(f"{t}: n={len(y):,} cases={int(sum(y)):,} | age beta={b0[1]:+.4f} | "
          f"age^2 beta={b1[2]:+.6f} | LRT chi2={lr:.2f} p={p:.3e}")
