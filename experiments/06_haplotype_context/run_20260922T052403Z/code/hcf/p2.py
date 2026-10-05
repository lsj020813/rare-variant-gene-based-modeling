
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import hashlib
import json
import os
import subprocess
import sys
import time

import numpy as np
import pandas as pd
from scipy import linalg as slinalg

from hcf import align as halign
from hcf import fastpath
from hcf import phase as hphase
from hcf import state as hstate
from hcf import tiles as htiles
from hcf import vcfio

BCFTOOLS = "bcftools"
ORIG = _config_path("${PROJECT_ROOT}/work/ref/orig_index/chr%s.vcf.gz")
RUN_ID = "HCF-20260922-v1/run_20260922T052403Z"
PROTOCOL_VERSION = "HCF-20260922-v1"

COV_COLS = ["age", "sex_male", "CT", "NC", "PC1", "PC2", "PC3", "PC4", "PC5"]
ALPHAS = (0.0, 0.5, 1.0)
LAMBDAS = (1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0)
INNER_FOLDS = 5
SEED = 20260922
MARKERS, STEP, MAX_GAP = 16, 8, 20000
NUISANCE_FLANK_BP = 500000
MIN_DISTINCT_A, MAX_STATES = 50, 32

MAX_FLANK_CANDIDATES = 400
MAX_NUISANCE_COLS = 256
NUIS_LD_R2 = 0.10
MAX_TILE_MARKERS = 900
MAX_DESIGN_COLS = 9000
TILE_PRUNE_R2_LADDER = (0.95, 0.90, 0.80, 0.60, 0.40, 0.20)
ENET_MAX_ITER = 3000
ENET_TOL = 1e-4
DF_PROBES = 20

ARMS = ("B_COV", "B_DS", "B_IND", "B_PAIR", "H_EXACT", "H_CLUSTER", "H_PHASE")
BASELINE_ARMS = ("B_DS", "B_IND", "B_PAIR")
CHALLENGER_ARMS = ("H_EXACT", "H_CLUSTER", "H_PHASE")
SIMPLICITY = {"B_DS": 1, "B_IND": 2, "B_PAIR": 3,
              "H_EXACT": 1, "H_CLUSTER": 2, "H_PHASE": 3}
ARM_BLOCKS = {
    "B_COV": (),
    "B_DS": ("NUIS", "DS"),
    "B_IND": ("NUIS", "DS", "IND"),
    "B_PAIR": ("NUIS", "DS", "IND", "PAIR"),
    "H_EXACT": ("NUIS", "DS", "IND", "ZEX"),
    "H_CLUSTER": ("NUIS", "DS", "IND", "ZCL"),
    "H_PHASE": ("NUIS", "DS", "IND", "PAIR", "Q"),
}

def sha256_text(s):
    return hashlib.sha256(s.encode("utf-8")).hexdigest()

def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def read_gt(chrom, region, samples_file, n_samples, keep):
    return vcfio.read_region(ORIG % chrom, region, n_samples=n_samples,
                             samples_file=samples_file, keep=keep)

def read_ds(chrom, region, samples_file, n_samples, keep=None, regions_file=None):
    cmd = [BCFTOOLS, "query"]
    if regions_file:
        cmd += ["-R", regions_file]
    else:
        cmd += ["-r", region]
    cmd += ["-S", samples_file, "-f", "%POS\t%REF\t%ALT[\t%DS]\n", ORIG % chrom]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=1 << 22)
    pos_l, ref_l, alt_l, rows = [], [], [], []
    try:
        for raw in proc.stdout:
            i1 = raw.index(b"\t")
            i2 = raw.index(b"\t", i1 + 1)
            i3 = raw.index(b"\t", i2 + 1)
            pos = int(raw[:i1])
            ref = raw[i1 + 1:i2].decode("ascii")
            alt = raw[i2 + 1:i3].decode("ascii")
            if keep is not None and (pos, ref, alt) not in keep:
                continue
            v = np.fromstring(raw[i3 + 1:].decode("ascii"), dtype=np.float64, sep="\t")
            if v.size != n_samples:
                raise vcfio.GtReadError(
                    "DS field count %d != n_samples %d at pos %d" % (v.size, n_samples, pos))
            pos_l.append(pos)
            ref_l.append(ref)
            alt_l.append(alt)
            rows.append(v.astype(np.float32))
    finally:
        try:
            proc.stdout.close()
        except Exception:
            pass
        proc.wait()
    ds = np.stack(rows, axis=0) if rows else np.zeros((0, n_samples), dtype=np.float32)
    return {"pos": np.asarray(pos_l, dtype=np.int64), "ref": ref_l, "alt": alt_l, "ds": ds}

def corr_matrix(X, idx):
    A = np.asarray(X[idx], dtype=np.float32)
    mu = A.mean(axis=0)
    sd = A.std(axis=0)
    sd[sd <= 0] = np.inf
    Z = (A - mu) / sd
    n = Z.shape[0]
    return (Z.T @ Z) / float(n)

def greedy_ld_prune(R, r2_max, force_keep=(), max_keep=None, order=None):
    m = R.shape[0]
    force = set(int(i) for i in force_keep)
    kept = []
    seq = list(order) if order is not None else list(range(m))
    for i in sorted(force):
        kept.append(i)
    for i in seq:
        if i in force:
            continue
        if max_keep is not None and len(kept) >= max_keep:
            break
        ok = True
        for j in kept:
            if R[i, j] * R[i, j] >= r2_max:
                ok = False
                break
        if ok:
            kept.append(i)
    return sorted(kept)

def encode_subwindow_keys(hap, subwins, verify=True):
    out = []
    for wi, (a, b) in enumerate(subwins):
        keys = fastpath.encode_keys(hap[:, :, a:b])
        if verify and wi == 0:
            fastpath.assert_equivalent_encoding(hap[:, :, a:b], keys)
        u, inv = np.unique(keys, return_inverse=True)
        uo = np.asarray([str(x) for x in u], dtype=object)
        out.append((uo, inv.reshape(keys.shape).astype(np.int32)))
        del keys
    return out

def transform_counts_codes(dictionary, uo, codes):
    cols = list(dictionary.columns)
    index = dict((c, j) for j, c in enumerate(cols))
    lut = np.array([index[dictionary.map_key(u)] for u in uo], dtype=np.int64)
    n = codes.shape[0]
    Z = np.zeros((n, len(cols)), dtype=np.float64)
    rows = np.arange(n)
    for cpy in range(2):
        np.add.at(Z, (rows, lut[codes[:, cpy]]), 1.0)
    if not np.all(Z.sum(axis=1) == 2.0):
        raise hstate.StateEncodingError("state count rows must sum to 2")
    return Z, cols

def state_block(keys_by_sw, train_idx, verify_first=False):
    mats, names, dicts = [], [], []
    for wi, (uo, codes) in enumerate(keys_by_sw):
        keys_obj = uo[codes]
        d = hstate.StateDictionary.fit(
            keys_obj[train_idx], min_distinct_people=MIN_DISTINCT_A, max_states=MAX_STATES,
            retain_other=True)
        if not d.states:
            dicts.append(d)
            continue
        Z, cols = transform_counts_codes(d, uo, codes)
        if verify_first and wi == 0:
            Zref, _c = fastpath.transform_counts(d, keys_obj)
            if not np.allclose(Z, Zref):
                raise AssertionError("code-path state counts differ from fastpath.transform_counts")
            fastpath.assert_equivalent_transform(d, keys_obj, Z)
        drop = cols.index(d.reference_state) if d.reference_state in cols else 0
        keepc = [j for j in range(len(cols)) if j != drop]
        if not keepc:
            dicts.append(d)
            continue
        mats.append(Z[:, keepc].astype(np.float32))
        names.extend(["w%03d|%s" % (wi, cols[j]) for j in keepc])
        dicts.append(d)
    if not mats:
        return np.zeros((keys_by_sw[0][1].shape[0], 0), dtype=np.float32), [], dicts
    return np.concatenate(mats, axis=1), names, dicts

def cluster_block(bca_train_map, n_people):
    mats, names = [], []
    for wi in range(bca_train_map.shape[0]):
        a = bca_train_map[wi].astype(np.int64)
        if a.size != 2 * n_people:
            raise ValueError("hapla bca haplotype count %d != 2*n %d" % (a.size, 2 * n_people))
        k = int(a.max()) + 1
        if k < 2:
            continue
        onehot = np.zeros((2 * n_people, k), dtype=np.float32)
        onehot[np.arange(2 * n_people), a] = 1.0
        Z = onehot[0::2] + onehot[1::2]
        ref = int(np.argmax(Z.sum(axis=0)))
        keepc = [j for j in range(k) if j != ref]
        mats.append(Z[:, keepc])
        names.extend(["c%03d|k%d" % (wi, j) for j in keepc])
    if not mats:
        return np.zeros((n_people, 0), dtype=np.float32), []
    return np.concatenate(mats, axis=1), names

def blocked_gram(X, rows, block=4096):
    p = X.shape[1]
    G = np.zeros((p, p), dtype=np.float64)
    for s in range(0, len(rows), block):
        idx = rows[s:s + block]
        Xb = np.asarray(X[idx], dtype=np.float64)
        G += Xb.T @ Xb
    return G

def enet_gram_solve(G, b, yvec, n, alpha, lam, w0=None):
    from sklearn.linear_model import _cd_fast
    p = G.shape[0]
    w = np.zeros(p, dtype=np.float64) if w0 is None else np.asarray(w0, dtype=np.float64).copy()
    l1 = float(alpha) * float(lam) * float(n)
    l2 = (1.0 - float(alpha)) * float(lam) * float(n)
    Q = np.asarray(G, dtype=np.float64, order="C")
    q = np.asarray(b, dtype=np.float64)
    rng = np.random.RandomState(SEED)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = _cd_fast.enet_coordinate_descent_gram(
            w, l1, l2, Q, q, np.asarray(yvec, dtype=np.float64),
            ENET_MAX_ITER, ENET_TOL, rng, 0, 0)
    gap, tol_reached = float(res[1]), float(res[2])
    converged = bool(gap <= tol_reached)
    return np.asarray(res[0], dtype=np.float64), int(res[3]), converged

def ridge_gram_solve(G, b, n, lam, want_df=False):
    p = G.shape[0]
    A = np.array(G, dtype=np.float64, order="C", copy=True)
    A.flat[::p + 1] += float(n) * float(lam)
    try:
        cf = slinalg.cho_factor(A, lower=True, check_finite=False)
    except slinalg.LinAlgError:
        A.flat[::p + 1] += 1e-8 * max(1.0, np.trace(G) / max(p, 1))
        cf = slinalg.cho_factor(A, lower=True, check_finite=False)
    beta = slinalg.cho_solve(cf, b, check_finite=False)
    df = float("nan")
    if want_df and p:
        rng = np.random.RandomState(SEED)
        z = rng.randint(0, 2, size=(p, DF_PROBES)).astype(np.float64) * 2.0 - 1.0
        Gz = G @ z
        S = slinalg.cho_solve(cf, Gz, check_finite=False)
        df = float(np.mean(np.einsum("ij,ij->j", z, S)))
    return beta, df

class FoldFit(object):

    def __init__(self, U, ncov, train_rows):
        self.U = U
        self.X = U[:, 1 + ncov:U.shape[1] - 1]
        self.y = U[:, -1]
        self.cov = U[:, 1:1 + ncov]
        self.train = np.asarray(train_rows)
        self.n = int(len(self.train))
        self.M = blocked_gram(U, self.train)
        self.ncov = ncov
        self.c_idx = np.arange(0, 1 + ncov)
        self.p = self.X.shape[1]
        Acc = self.M[np.ix_(self.c_idx, self.c_idx)]
        self.Acc_inv = np.linalg.pinv(Acc)
        self.Acy = self.M[self.c_idx, -1]
        self.gamma_y = self.Acc_inv @ self.Acy
        self.yy = float(self.M[-1, -1])
        self.sse_cov = max(self.yy - float(self.Acy @ self.gamma_y), 0.0)
        diag = np.diag(self.M)[1 + ncov:1 + ncov + self.p]
        mean = self.M[0, 1 + ncov:1 + ncov + self.p] / self.n
        var = diag / self.n - mean * mean
        var[var < 0] = 0.0
        self.sd = np.sqrt(var)
        self.y_res_train = None

    def y_resid(self):
        if self.y_res_train is None:
            C = np.empty((self.n, 1 + self.ncov), dtype=np.float64)
            C[:, 0] = 1.0
            C[:, 1:] = self.cov[self.train]
            self.y_res_train = self.y[self.train].astype(np.float64) - C @ self.gamma_y
        return self.y_res_train

    def arm_system(self, cols):
        cols = np.asarray(cols, dtype=np.int64)
        if cols.size == 0:
            return None
        keep = cols[self.sd[cols] > 0]
        if keep.size == 0:
            return None
        gi = keep + 1 + self.ncov
        Axx = self.M[np.ix_(gi, gi)]
        Acx = self.M[np.ix_(self.c_idx, gi)]
        Axy = self.M[gi, -1]
        W = self.Acc_inv @ Acx
        G = Axx - Acx.T @ W
        b = Axy - W.T @ self.Acy
        s = self.sd[keep]
        D = 1.0 / s
        G = G * D[:, None] * D[None, :]
        b = b * D
        return G, b, keep, s

    def solve(self, cols, alpha, lam, want_df=False, w0=None):
        sysm = self.arm_system(cols)
        if sysm is None:
            return {"beta_raw": np.zeros(0), "keep": np.zeros(0, dtype=np.int64),
                    "gamma": self.gamma_y, "df": 0.0, "nnz": 0, "beta_std": np.zeros(0),
                    "converged": True}
        G, b, keep, s = sysm
        converged = True
        if alpha == 0.0:
            beta_std, df = ridge_gram_solve(G, b, self.n, lam, want_df=want_df)
            nnz = int(np.sum(np.abs(beta_std) > 0))
        else:
            beta_std, _ni, converged = enet_gram_solve(G, b, self.y_resid(), self.n,
                                                       alpha, lam, w0=w0)
            nnz = int(np.sum(beta_std != 0))
            df = float(nnz)
        beta_raw = beta_std / s
        gamma = self.gamma_y - (self.Acc_inv @ (self.M[np.ix_(self.c_idx, keep + 1 + self.ncov)]
                                                @ beta_raw))
        return {"beta_raw": beta_raw, "keep": keep, "gamma": gamma, "df": df,
                "nnz": nnz, "beta_std": beta_std, "converged": bool(converged)}

    def predict(self, fit, rows):
        rows = np.asarray(rows)
        C = np.empty((len(rows), 1 + self.ncov), dtype=np.float64)
        C[:, 0] = 1.0
        C[:, 1:] = self.cov[rows]
        yh = C @ fit["gamma"]
        if fit["keep"].size:
            yh = yh + np.asarray(self.X[np.ix_(rows, fit["keep"])], dtype=np.float64) @ fit["beta_raw"]
        return yh

def one_se_select(cv_mse, cv_se_at_min, cands):
    j = int(np.argmin(cv_mse))
    thr = cv_mse[j] + cv_se_at_min
    ok = [i for i in range(len(cands)) if cv_mse[i] <= thr]
    ok.sort(key=lambda i: (-cands[i][1], -cands[i][0]))
    return ok[0], j

def paired_gain(y, yh_base, yh_chal, s2):
    d = ((y - yh_base) ** 2 - (y - yh_chal) ** 2) / s2
    n = d.size
    m = float(d.mean())
    sd = float(d.std(ddof=1)) if n > 1 else float("nan")
    se = sd / np.sqrt(n) if n > 1 else float("nan")
    z = m / se if se and np.isfinite(se) and se > 0 else float("nan")
    from scipy import stats as sstats
    p = float(sstats.norm.sf(z)) if np.isfinite(z) else float("nan")
    return m, se, z, p
