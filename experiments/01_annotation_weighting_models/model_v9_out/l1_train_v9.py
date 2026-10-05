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
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys
import time

OWNER = Path(_config_path("${PROJECT_ROOT}"))
REF = OWNER / "work/ref15"
PREREG = OWNER / "work/run_band15/traits_v8_prereg.tsv"
FOLDS = {0: {"7", "13", "16", "19", "21"},
         1: {"2", "8", "9", "10", "15"}, 2: {"1", "6", "14"},
         3: {"3", "12", "17", "22"}, 4: {"4", "5", "11", "18", "20"}}
LEGACY = ("tchl", "htn", "dm", "lip")
S1 = ["cons", "epi_active", "epi_repr", "epi_trans", "tf", "linsight",
      "gpn_msa", "dist_tss"]
S2 = ["re2g_max", "gh_elem_score", "gh_link_score"]
S3 = ["is_cage_prom", "in_body", "in_tss3kb", "in_re2g", "t1_na", "is_indel"]
DIAG = ["maf", "r2", "avg_cs"]
STAAR_ANNOTATIONS = ["cadd", "cons", "epi_active", "epi_repr", "epi_trans",
                     "tf", "linsight", "gpn_msa"]
DEFAULT_KAPPA_GRID = [0.0, 0.25, 0.5, 0.75, 1.0]
DEFAULT_COND_CHECK_GENES = ["LDLR", "HMGCR", "PCSK9", "APOB"]
CL = 2.0
EPS_F = 0.05
NKNOT = 6
NUM_EPS = 1e-8
ZERO_SD = 1e-12
T0 = time.time()

def require(condition, message):
    if not condition:
        raise RuntimeError("GATE FAIL: " + message)

def log(message):
    print(f"[{time.time() - T0:7.0f}s] {message}", flush=True)

def sha_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True,
                                    allow_nan=False).encode()).hexdigest()

def owned_path(path):
    p = Path(path).resolve()
    try:
        p.relative_to(OWNER)
    except ValueError:
        require(False, "write path must be under owner directory")
    parent = p if p.exists() else p.parent
    while not parent.exists():
        parent = parent.parent
    require(parent.stat().st_uid == os.getuid(), "write ancestor owner mismatch")
    return p

def atomic_json(path, value):
    path = owned_path(path)
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "w") as fh:
        json.dump(value, fh, ensure_ascii=False, indent=2, allow_nan=False)
        fh.write("\n")
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)

def atomic_tsv(path, fieldnames, rows):
    path = owned_path(path)
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames,
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)

def atomic_npz(path, **arrays):
    path = owned_path(path)
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "wb") as fh:
        np.savez(fh, **arrays)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)

def gene_key(value):
    return value.strip().split(".")[0]

def norm_iid(value):
    value = value.strip()
    parts = value.split("_")
    return parts[0] if len(parts) == 2 and parts[0] == parts[1] else value

def read_traits(path):
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = set(reader.fieldnames or [])
        require(len(fields) == len(reader.fieldnames or []), "duplicate trait TSV headers")
        require({"trait", "type", "transform"} <= fields, "trait TSV schema")
        wide = {"CT_col", "NC_col", "AS_col", "table"} <= fields
        require(wide or {"file", "col"} <= fields, "trait TSV source columns")
        rows = []
        for record in reader:
            require(None not in record and all(v is not None for v in record.values()), "trait TSV row width")
            rows.append({k: v.strip() for k, v in record.items()})
    require(bool(rows), "empty trait list")
    names = [r["trait"] for r in rows]
    require(len(names) == len(set(names)), "duplicate trait identifiers")
    for row in rows:
        require(re.fullmatch(r"[a-z][a-z0-9_]*", row["trait"]) is not None,
                "trait identifiers must be safe lowercase names")
        require(row["type"] in ("quant", "binary"), "unknown trait type")
        require(row["trait"] != "mean", "mean is reserved for aggregate metrics")
        require(row["transform"] in ("rint", "rint(log)", "-", "none"),
                "unsupported transform; no silent transform substitution")
    return rows, "cohort_columns" if wide else "file_column"

def safety_args(parser):
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--memory-gb", type=float, default=120)
    parser.add_argument("--tmpdir", default=str(OWNER / "work/tmp"))
    parser.add_argument("--device", choices=["cpu", "cuda"], default="cuda")

def resource_guard(args):
    sys.dont_write_bytecode = True
    import fcntl
    import resource
    import subprocess
    require(1 <= args.threads <= 8, "threads must be 1..8")
    require(0 < args.memory_gb <= 120, "memory cap must be <=120 GiB")
    tmp = owned_path(args.tmpdir)
    require(tmp.is_dir(), "create owner TMPDIR before launching")
    for name in ("TMPDIR", "TMP", "TEMP"):
        os.environ[name] = str(tmp)
    os.environ["XDG_CACHE_HOME"] = str(tmp / "cache")
    os.environ["CUDA_CACHE_PATH"] = str(tmp / "cuda-cache")
    os.environ["TORCH_HOME"] = str(tmp / "torch")
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                 "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[name] = str(args.threads)
    os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
    cap = int(args.memory_gb * 1024 ** 3)
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    cap = min([cap] + [x for x in (soft, hard) if x != resource.RLIM_INFINITY])
    resource.setrlimit(resource.RLIMIT_AS, (cap, cap))
    meminfo = dict(line.split(":", 1) for line in Path("/proc/meminfo").read_text().splitlines())
    available = int(meminfo["MemAvailable"].split()[0]) * 1024
    require(available >= 80 * 1024 ** 3, "available RAM below policy 80 GiB")
    heavy = ("deepsea", "deepripe", "saige", "step2", "deeprvat", "sparse_geno",
             "l1_train", "build_resid", "plink2")
    ancestors = {os.getpid()}
    parent = os.getppid()
    while parent > 0 and parent not in ancestors:
        ancestors.add(parent)
        try:
            fields = Path(f"/proc/{parent}/stat").read_text().rsplit(")", 1)[1].split()
            parent = int(fields[1])
        except FileNotFoundError:
            break
    def snapshot():
        procs = {}
        for p in Path("/proc").glob("[0-9]*"):
            try:
                pid = int(p.name)
                uid = p.stat().st_uid
                stat = (p / "stat").read_text().rsplit(")", 1)[1].split()
                ticks = int(stat[11]) + int(stat[12])
                rss = int(stat[21]) * os.sysconf("SC_PAGE_SIZE")
                procs[pid] = (uid, ticks, rss, stat[19])
                if uid == os.getuid() and pid not in ancestors:
                    tokens = (p / "cmdline").read_bytes().decode(errors="replace").lower().split("\0")
                    tokens = [t for t in tokens if t]
                    viewer = {"tail", "cat", "grep", "rg", "less", "more", "head", "sed", "awk", "wc", "ls", "vim", "nano", "stat", "pgrep", "tmux"}
                    if tokens and os.path.basename(tokens[0]) not in viewer and stat[0] not in ("Z", "z", "X", "x"):
                        basenames = [os.path.basename(t.split("=")[-1]) for t in tokens]
                        require(not any(b.startswith(heavy) for b in basenames), "another owner heavy stage is active")
            except (FileNotFoundError, ProcessLookupError):
                continue
            except PermissionError:
                require(False, "cannot observe process usage; use an approved visible guard")
        return procs
    start = time.monotonic()
    first = snapshot()
    time.sleep(1.0)
    second = snapshot()
    elapsed = time.monotonic() - start
    users = {}
    for pid, (uid, ticks, rss, born) in second.items():
        if uid == os.getuid():
            continue
        usage = users.setdefault(uid, [0.0, 0])
        old = first.get(pid)
        delta = ticks - old[1] if old and old[3] == born else 0
        usage[0] += max(0, delta) / os.sysconf("SC_CLK_TCK") / elapsed * 100
        usage[1] += rss
    require(all(v[0] <= 400 and v[1] <= 32 * 1024 ** 3 for v in users.values()),
            "another user's CPU/RAM exceeds policy")
    require(sum(v[0] for v in users.values()) <= 800 and
            sum(v[1] for v in users.values()) <= 64 * 1024 ** 3,
            "aggregate other-user CPU/RAM exceeds policy")
    if args.device == "cuda":
        reply = subprocess.run(["nvidia-smi", "--query-compute-apps=pid",
                                "--format=csv,noheader,nounits"],
                               capture_output=True, text=True, check=True)
        require(not reply.stdout.strip(), "GPU already has a compute process")
    lock = open(tmp / "l1_v8_heavy.lock", "a+")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        lock.close()
        require(False, "another v8/v9 trainer or builder holds the heavy-stage lock")
    owned_path(args.out).mkdir(parents=True, exist_ok=True)
    return lock, dict(threads=args.threads, memory_cap_bytes=cap, tmpdir=str(tmp),
                      available_ram_bytes=available, other_users=len(users),
                      guard="v9 stdlib usage guard; external stages detected by command names")

def load_libraries(threads):
    global np, sp, torch, minimize, null_space, SplineTransformer, norm, rankdata
    import numpy as np
    import scipy.sparse as sp
    import torch
    from scipy.optimize import minimize
    from scipy.linalg import null_space
    from scipy.stats import norm, rankdata
    from sklearn.preprocessing import SplineTransformer
    torch.set_num_threads(threads)
    torch.set_num_interop_threads(1)
    torch.backends.cuda.matmul.allow_tf32 = False

class Data:
    def __init__(self, args):
        self.args = args
        self.root = Path(args.root)
        self.fingerprints = {}
        path = self.root / "annot/cache/fm_all.npz"
        self.fingerprint(path)
        with np.load(path, allow_pickle=True) as c:
            self.X = c["X"]
            self.cols = list(c["cols"].astype(str))
            self.keys = c["key37"].astype(str)
            self.genes = c["gene"].astype(str)
            self.chrs = c["chr"].astype(int).astype(str)
        summaries = sorted((self.root / "groupfiles_bwg").glob("chr*.summary.json"))
        require(bool(summaries), "group summary files absent")
        expected = 0
        for path in summaries:
            with open(path) as fh:
                expected += int(json.load(fh)["pairs"])
        require(self.X.shape == (expected, 33), "cache shape vs group summaries")
        require(len(set(self.cols)) == len(self.cols), "duplicate cache columns")
        require(set(S1 + S2 + S3 + DIAG + ["n_genes"]) <= set(self.cols), "cache columns")
        require(all(len(x) == expected for x in (self.keys, self.genes, self.chrs)),
                "cache annotation lengths")
        if args.smoke:
            sel = self.chrs == args.smoke
            self.X = self.X[sel]
            self.keys = self.keys[sel]
            self.genes = self.genes[sel]
            self.chrs = self.chrs[sel]
        self.phi_continuous = list(S1)
        self.cadd_join = None
        if args.v9_enabled:
            self.load_cadd(args.cadd_side)
        self.use_chr = sorted(set(self.chrs), key=int)
        require(self.use_chr == ([args.smoke] if args.smoke else [str(i) for i in range(1, 23)]),
                "expected chromosomes missing")
        self.gl = sorted(set(self.genes))
        canonical = [gene_key(g) for g in self.gl]
        if args.v9_enabled:
            require(len(canonical) == len(set(canonical)),
                    "version-stripped cache gene IDs are not unique")
        self.gix = {g: i for i, g in enumerate(self.gl)}
        self.gix_key = {gene_key(g): i for i, g in enumerate(self.gl)}
        self.gid = np.array([self.gix[g] for g in self.genes], dtype=np.int32)
        self.ng = len(self.gl)
        gcs = {}
        for g, c in zip(self.genes, self.chrs):
            gcs.setdefault(g, set()).add(c)
        require(all(len(v) == 1 for v in gcs.values()), "multichromosome gene identifiers")
        self.gchr = {g: next(iter(cs)) for g, cs in gcs.items()}
        order = np.argsort(self.gid, kind="stable")
        cuts = np.searchsorted(self.gid[order], np.arange(self.ng + 1))
        self.pairs_of = [order[cuts[g]:cuts[g + 1]] for g in range(self.ng)]
        require(sum(map(len, self.pairs_of)) == len(self.keys), "gene pair partition")
        n_genes = self.col("n_genes")
        require(np.isfinite(n_genes).all() and (n_genes >= 1).all(), "invalid n_genes")
        self.pw = (1 / n_genes).astype(np.float32)
        self.dsc = {}
        self.cond = None
        if not args.export_basis_only:
            self.load_ds()
            if args.v9_enabled:
                self.load_conditioning(args.cond_cov_path)

    def fingerprint(self, path):
        st = path.stat()
        self.fingerprints[str(path)] = dict(size=st.st_size, mtime_ns=st.st_mtime_ns)

    def col(self, name):
        return self.X[:, self.cols.index(name)].astype(np.float64)

    def load_cadd(self, path):
        path = Path(path)
        self.fingerprint(path)
        side = {}
        rows = 0
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            require((reader.fieldnames or []) == ["variant_key", "cadd_phred"],
                    "cadd_side.tsv schema/order must be variant_key,cadd_phred")
            for record in reader:
                rows += 1
                key = record["variant_key"].strip()
                require(bool(key), "empty CADD variant_key")
                require(key not in side, "duplicate CADD variant_key")
                value_text = record["cadd_phred"].strip()
                value = np.nan if value_text.lower() in ("", ".", "na", "nan") else float(value_text)
                require(np.isnan(value) or (math.isfinite(value) and value >= 0), "invalid CADD PHRED")
                side[key] = value
        require(rows > 0, "empty cadd_side.tsv")
        key_hit = np.asarray([k in side for k in self.keys], dtype=bool)
        require(key_hit.any(), "CADD join matched zero variant keys")
        cadd = np.array([side.get(k, np.nan) for k in self.keys], dtype=np.float64)
        matched = np.isfinite(cadd)
        require(matched.any(), "CADD join matched zero cache rows")
        t1_na = self.col("t1_na").astype(bool)
        require(np.array_equal(~matched, t1_na),
                "CADD missingness must equal the cache t1_na convention")
        self.X = np.column_stack([self.X, cadd]).astype(self.X.dtype, copy=False)
        self.cols.append("cadd")
        self.phi_continuous.append("cadd")
        self.cadd_join = dict(side_unique_keys=len(side), cache_pair_rows=len(self.keys),
                              key_matched_pair_rows=int(key_hit.sum()),
                              matched_pair_rows=int(matched.sum()), missing_pair_rows=int((~matched).sum()),
                              side_duplicate_keys=0, join_key="variant_key",
                              missing_rule="missing iff t1_na=1; centered-zero continuous value plus t1_na indicator")

    def load_conditioning(self, path):
        path = Path(path)
        self.fingerprint(path)
        with np.load(path, allow_pickle=False) as z:
            require(set(z.files) == {"ids", "X", "snp_list"},
                    "cond_covariates.npz members must be exactly ids,X,snp_list")
            ids = z["ids"].astype(str)
            X = z["X"]
            snp_list = z["snp_list"].astype(str)
        require(X.ndim == 2 and ids.ndim == 1 and snp_list.ndim == 1,
                "conditioning array dimensions")
        require(X.shape == (len(ids), len(snp_list)), "conditioning X shape")
        require(np.issubdtype(X.dtype, np.number) and np.isfinite(X).all(),
                "conditioning X must be finite numeric")
        normalized = [norm_iid(x) for x in ids]
        require(len(normalized) == len(set(normalized)), "duplicate normalized conditioning ID")
        require(len(snp_list) == len(set(snp_list)), "duplicate conditioning snp_list key")
        wanted = [norm_iid(x) for x in self.samples]
        require(len(wanted) == len(set(wanted)), "duplicate normalized DS ID")
        loc = {key: i for i, key in enumerate(normalized)}
        require(len(loc) == len(wanted) and set(loc) == set(wanted),
                "conditioning/DS ID sets differ")
        order = np.asarray([loc[key] for key in wanted], dtype=np.int64)
        X = np.ascontiguousarray(X[order], dtype=np.float64)
        if X.shape[1]:
            xtx = X.T @ X
            require(np.linalg.matrix_rank(xtx) == X.shape[1],
                    "conditioning X'X is singular")
            xtx_inv = np.linalg.inv(xtx)
            require(np.isfinite(xtx_inv).all(), "nonfinite conditioning inverse")
        else:
            xtx_inv = np.empty((0, 0), dtype=np.float64)
        self.cond = dict(X=X, xtx_inv=xtx_inv, n_covariates=X.shape[1],
                         snp_list_sha256=digest(snp_list.tolist()),
                         id_join=dict(ds_unique=len(wanted), cov_unique=len(normalized),
                                      ds_duplicates=0, cov_duplicates=0, unmatched=0,
                                      reordered=bool(np.any(order != np.arange(len(order)))),
                                      join_key="normalized ids"),
                         inverse_calculations=1)
        log(f"conditional covariates loaded: N={len(ids)}, C={X.shape[1]}; XTX inverse cached once")

    def load_ds(self):
        self.vloc = np.empty(len(self.keys), np.int32)
        for chrom in self.use_chr:
            path = self.root / f"annot/ds/chr{chrom}.ds.npz"
            self.fingerprint(path)
            with np.load(path, allow_pickle=True) as z:
                samples = z["samples"].astype(str).tolist()
                require(len(samples) == N_SAMPLES or self.args.synthetic_smoke,
                        "sample count differs from v72 cohort")
                require(len(set(map(norm_iid, samples))) == len(samples), "normalized sample key duplicates")
                if not self.dsc:
                    self.samples = samples
                    self.ns = len(samples)
                else:
                    require(samples == self.samples, "chromosome sample order mismatch")
                indptr = z["indptr"]
                indices = z["indices"]
                values = z["data"]
                variant_keys = z["keys"].astype(str)
                require(0 <= int(indptr[-1]) < 2 ** 31, "chromosome nnz exceeds int32")
                require(len(variant_keys) < 2 ** 31, "chromosome rows exceed int32")
                require(values.dtype == np.float32, "DS data must be prebuilt float32")
                require(len(indices) == len(values) == int(indptr[-1]), "DS nnz mismatch")
                require(len(indptr) == len(variant_keys) + 1, "DS indptr length")
                require(indptr[0] == 0 and np.all(indptr[1:] >= indptr[:-1]), "DS indptr monotonicity")
                require(indices.min() >= 0 and indices.max() < self.ns, "DS column bounds")
                matrix = sp.csr_matrix((values, indices.astype(np.int32, copy=False),
                                       indptr.astype(np.int32, copy=False)),
                                      shape=(len(variant_keys), self.ns), copy=False)
                require(matrix.indices.dtype == matrix.indptr.dtype == np.int32, "CSR index dtype")
                row_of = {key: i for i, key in enumerate(variant_keys)}
                require(len(row_of) == len(variant_keys), "DS variant key duplicates")
                pair_idx = np.flatnonzero(self.chrs == chrom)
                require(all(k in row_of for k in self.keys[pair_idx]), "annotation to DS join misses")
                self.vloc[pair_idx] = [row_of[k] for k in self.keys[pair_idx]]
                self.dsc[chrom] = matrix
                del row_of, values, indices, indptr, variant_keys
            log(f"DS chromosome {chrom} loaded; int32 indices verified")
        for g, pp in enumerate(self.pairs_of):
            require(np.all(self.chrs[pp] == self.gchr[self.gl[g]]), "gene/DS chromosome mismatch")
            require(len(set(self.keys[pp])) == len(pp), "duplicate variant/gene annotation key")

    def variant_betas(self, matrix, rows):
        if self.cond is None or self.cond["n_covariates"] == 0:
            return np.empty((len(rows), 0), dtype=np.float64)
        return (matrix[rows] @ self.cond["X"]) @ self.cond["xtx_inv"]

    def burden(self, phi, gi, conditioned=True):
        result = np.empty((len(gi), self.ns), dtype=phi.dtype)
        for row, g in enumerate(gi):
            pp = self.pairs_of[g]
            matrix = self.dsc[self.gchr[self.gl[g]]]
            rows = self.vloc[pp]
            weights = phi[pp]
            raw = matrix[rows].T @ weights
            if conditioned and self.cond is not None and self.cond["n_covariates"]:
                beta = self.variant_betas(matrix, rows)
                raw = raw - self.cond["X"] @ (weights.astype(np.float64) @ beta)
            result[row] = raw
        return result

    def grad_phi(self, gS, gi, result, conditioned=True):
        for row, g in enumerate(gi):
            pp = self.pairs_of[g]
            matrix = self.dsc[self.gchr[self.gl[g]]]
            rows = self.vloc[pp]
            value = matrix[rows] @ gS[row]
            if conditioned and self.cond is not None and self.cond["n_covariates"]:
                beta = self.variant_betas(matrix, rows)
                value = value - beta @ (self.cond["X"].T @ gS[row].astype(np.float64))
            result[pp] = value

class Residuals:
    def __init__(self, args, data, rows):
        self.names = [row["trait"] for row in rows]
        self.full = np.full((len(rows), data.ns), np.nan, np.float32)
        self.hashes = {}
        self.kind = {}
        self.summary = {}
        index = {norm_iid(s): i for i, s in enumerate(data.samples)}
        if args.resid == "cov":
            root = Path(args.resid_dir)
            summary_path = root / "build_resid_v8.summary.json"
            with open(summary_path) as fh:
                self.summary = json.load(fh)
            with open(root / "build_resid_v8.done") as fh:
                marker = json.load(fh)
            require(marker["summary_sha256"] == sha_file(summary_path), "residual completion marker")
            require(self.summary["traits_sha256"] == sha_file(args.traits), "residual/prereg hash mismatch")
            require(self.summary["sample_order_sha256"] == digest(data.samples), "builder/genotype sample order mismatch")
        for it, row in enumerate(rows):
            trait = row["trait"]
            if args.resid == "eta":
                require(trait in LEGACY, "eta supports only the four legacy traits")
                path = Path(args.offset_dir) / f"{trait}.eta.tsv"
                self.kind[trait] = "eta_y_minus_eta" if row["type"] == "quant" else "eta_y_minus_mu_working"
            else:
                path = Path(args.resid_dir) / f"{trait}.resid.tsv"
                self.kind[trait] = "cov_ols_age_sex_cohort_PC1to10_then_RINT_supervision_only"
            self.hashes[trait] = sha_file(path)
            if args.resid == "cov":
                require(self.hashes[trait] == self.summary["traits"][trait]["residual_sha256"],
                        "residual content differs from builder manifest")
            seen = set()
            matched = 0
            with open(path, newline="") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                fields = reader.fieldnames or []
                require(bool(fields), "empty residual input")
                idcol = "IID" if args.resid == "cov" else fields[0]
                needed = {idcol, "resid"} if args.resid == "cov" else {idcol, "y", "eta", "mu"}
                require(needed <= set(fields), "residual schema mismatch")
                for record in reader:
                    key = norm_iid(record[idcol])
                    require(key not in seen, "duplicate normalized residual IID")
                    seen.add(key)
                    if key not in index:
                        continue
                    if args.resid == "cov":
                        value = float(record["resid"])
                    else:
                        offset = "eta" if row["type"] == "quant" else "mu"
                        value = float(record["y"]) - float(record[offset])
                    if math.isfinite(value):
                        self.full[it, index[key]] = value
                        matched += 1
            require(matched >= args.min_n, f"trait {trait}: insufficient matched residuals")
        if args.perm_r:
            perm = np.random.default_rng(args.perm_r).permutation(data.ns)
            self.full = self.full[:, perm]
        self.mask = np.isfinite(self.full)
        self.R = np.zeros_like(self.full, dtype=np.float32)
        self.n = self.mask.sum(axis=1)
        for it in range(len(rows)):
            observed = self.full[it, self.mask[it]].astype(np.float64)
            sd = observed.std(ddof=0)
            require(sd > NUM_EPS, "zero residual variance")
            self.R[it, self.mask[it]] = (observed - observed.mean()) / sd / np.sqrt(self.n[it])
        packed = np.packbits(self.mask.T, axis=1)
        groups = {}
        strata = [""] * data.ns
        if args.perm_strata:
            seen = set()
            with open(args.perm_strata, newline="") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                require({"IID", "stratum"} <= set(reader.fieldnames or []), "permutation strata schema")
                for record in reader:
                    key = norm_iid(record["IID"])
                    require(key not in seen, "duplicate permutation-stratum IID")
                    seen.add(key)
                    if key in index:
                        require(bool(record["stratum"]), "empty stratum")
                        strata[index[key]] = record["stratum"]
            require(all(strata), "strata missing for genotype samples")
            if args.perm_r:
                strata = [strata[i] for i in perm]
        for i, bits in enumerate(packed):
            groups.setdefault((bits.tobytes(), strata[i]), []).append(i)
        self.groups = [np.asarray(v, np.int32) for v in groups.values()]
        self.perm_meta = dict(kind="joint within exact observed-trait mask and optional stratum",
                              n_blocks=len(self.groups), singleton_samples=sum(len(v) == 1 for v in self.groups),
                              strata_sha256=sha_file(args.perm_strata) if args.perm_strata else None)
        movable = sum(len(v) for v in self.groups if len(v) > 1 and self.mask[:, v[0]].any())
        require(movable >= args.min_n, "too few exchangeable observed individuals")
        self.perm_meta["movable_observed_samples"] = movable

    def permutations(self, count, seed):
        rng = np.random.default_rng(seed)
        permutations = []
        for _ in range(count):
            perm = np.arange(self.R.shape[1], dtype=np.int32)
            for group in self.groups:
                perm[group] = rng.permutation(group)
            permutations.append(perm)
        return permutations

    def diagnostics(self):
        nt = len(self.names)
        corr = np.eye(nt)
        overlap = np.zeros((nt, nt), dtype=int)
        for i in range(nt):
            for j in range(i, nt):
                mask = self.mask[i] & self.mask[j]
                overlap[i, j] = overlap[j, i] = int(mask.sum())
                if mask.sum() < 2:
                    corr[i, j] = corr[j, i] = np.nan
                else:
                    x = self.full[i, mask].astype(np.float64)
                    y = self.full[j, mask].astype(np.float64)
                    value = np.corrcoef(x, y)[0, 1] if x.std() > 0 and y.std() > 0 else np.nan
                    corr[i, j] = corr[j, i] = value
        gram = self.R.astype(np.float64) @ self.R.astype(np.float64).T
        ev = np.linalg.eigvalsh((gram + gram.T) / 2)
        effective = float(ev.sum() ** 2 / np.square(ev).sum())
        cc = int(self.mask.all(axis=0).sum())
        return dict(traits=self.names, n={t: int(n) for t, n in zip(self.names, self.n)},
                    common_genotype_n=self.R.shape[1], complete_case_n=cc,
                    intersection_loss={t: int(n - cc) for t, n in zip(self.names, self.n)},
                    pairwise_correlation=[[float(x) if np.isfinite(x) else None for x in row] for row in corr],
                    pairwise_overlap=overlap.tolist(), supervision_gram_correlation=gram.tolist(),
                    eigenvalues=ev.tolist(), effective_traits=effective,
                    effective_formula="(sum eigenvalues)^2 / sum(eigenvalues^2), on PSD zero-filled normalized residual Gram",
                    effective_status="descriptive participation ratio; no loss weighting or multiple-testing correction",
                    R_dtype="float32", R_bytes=self.R.nbytes, resid_kind=self.kind,
                    supervision_only=True, permutation=self.perm_meta)

class BBJTruth:
    def __init__(self, args, trait_names):
        self.z_path = Path(args.bbj_gene_z)
        self.set_path = Path(args.bbj_gene_set)
        self.traits = list(trait_names)
        self.z = {}
        self.n_snps = {}
        with open(self.z_path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            require((reader.fieldnames or []) == ["gene_id", "trait", "z_cond", "n_snps"],
                    "bbj_gene_z.tsv schema/order must be gene_id,trait,z_cond,n_snps")
            for record in reader:
                key = (gene_key(record["gene_id"]), record["trait"].strip())
                require(key[0] and key[1], "empty BBJ gene/trait key")
                require(key not in self.z, "duplicate BBJ gene_id/trait key")
                value = float(record["z_cond"])
                require(math.isfinite(value), "nonfinite BBJ z_cond")
                n_text = record["n_snps"].strip()
                require(re.fullmatch(r"[0-9]+", n_text) is not None, "BBJ n_snps must be an integer")
                n_snps = int(n_text)
                require(n_snps > 0, "BBJ n_snps must be positive")
                self.z[key] = value
                self.n_snps[key] = n_snps
        require(bool(self.z), "empty bbj_gene_z.tsv")
        require({t for _, t in self.z} == set(self.traits),
                "BBJ z traits must exactly match the preregistered trait list")
        self.sig = {}
        with open(self.set_path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            require((reader.fieldnames or []) == ["gene_id", "trait", "is_bbj_sig"],
                    "bbj_gene_set.tsv schema/order must be gene_id,trait,is_bbj_sig")
            for record in reader:
                key = (gene_key(record["gene_id"]), record["trait"].strip())
                require(key[0] and key[1], "empty BBJ-set gene/trait key")
                require(key not in self.sig, "duplicate BBJ-set gene_id/trait key")
                require(record["is_bbj_sig"].strip() in ("0", "1"),
                        "is_bbj_sig must be 0 or 1")
                self.sig[key] = record["is_bbj_sig"].strip() == "1"
        require(set(self.sig) == set(self.z),
                "BBJ z/set gene_id,trait key sets must match exactly")
        self.target = {}
        self.z_standard = {}
        self.trait_stats = {}
        for trait in self.traits:
            keys = [key for key in self.z if key[1] == trait]
            values = np.asarray([abs(self.z[key]) for key in keys], dtype=np.float64)
            mu = float(values.mean())
            sigma = float(values.std(ddof=0))
            require(sigma > NUM_EPS, f"BBJ trait {trait}: zero |z_cond| SD")
            standardized = (values - mu) / sigma
            upper = float(np.quantile(standardized, args.target_clip_q))
            clipped = np.minimum(standardized, upper)
            for key, raw, value in zip(keys, standardized, clipped):
                self.z_standard[key] = float(raw)
                self.target[key] = float(value)
            self.trait_stats[trait] = dict(n=len(keys), abs_z_mean=mu, abs_z_sd=sigma,
                                           target_clip_q=args.target_clip_q,
                                           target_clip_upper=upper,
                                           n_clipped=int(np.sum(standardized > upper)))
        self.sigmoid_scale = args.bbj_prior_sigmoid_scale
        self.join_audit = dict(z_key="version-stripped gene_id + trait",
                               set_key="version-stripped gene_id + trait",
                               z_duplicate_keys=0, set_duplicate_keys=0,
                               key_set_symmetric_difference=0)

    def validate_against_data(self, data, active_traits, minimum):
        missing = sum((gene_key(g), trait) not in self.z
                      for g in data.gl for trait in active_traits)
        designated = {gene for (gene, trait), flag in self.sig.items()
                      if flag and trait in active_traits and gene in data.gix_key}
        require(len(designated) >= minimum,
                f"U17 needs at least {minimum} BBJ-designated genes with band variants; observed {len(designated)}")
        self.join_audit.update(cache_unique_gene_keys=data.ng, active_traits=len(active_traits),
                               missing_cache_gene_trait_keys=int(missing),
                               designated_band_genes=len(designated), min_bbj_genes=minimum)

    def arrays(self, data, gi, names):
        genes = [gene_key(data.gl[g]) for g in gi]
        tested = np.asarray([[(gene, trait) in self.target for trait in names]
                             for gene in genes], dtype=bool)
        targets = np.asarray([[self.target.get((gene, trait), np.nan) for trait in names]
                              for gene in genes], dtype=np.float64)
        sig = np.asarray([[self.sig.get((gene, trait), False) for trait in names]
                          for gene in genes], dtype=bool)
        zstd = np.asarray([[self.z_standard.get((gene, trait), np.nan) for trait in names]
                           for gene in genes], dtype=np.float64)
        scaled = self.sigmoid_scale * zstd
        per_trait_prior = np.exp(-np.logaddexp(0.0, -scaled))
        counts = tested.sum(axis=1)
        prior = np.divide(np.nansum(per_trait_prior, axis=1), counts,
                          out=np.full(len(genes), np.nan), where=counts > 0)
        observed_prior = prior[np.isfinite(prior)]
        require(np.isfinite(targets[tested]).all() and
                ((observed_prior > 0) & (observed_prior < 1)).all(),
                "invalid BBJ target/prior arrays")
        return targets, sig, prior, tested

    def manifest(self):
        return dict(target_definition="clip((abs(z_cond)-trait_mean)/trait_sd, upper trait quantile)",
                    target_clip="upper only", prior_definition="trait mean sigmoid(scale * standardized abs(z_cond))",
                    prior_sigmoid_scale=self.sigmoid_scale, trait_stats=self.trait_stats,
                    join_audit=self.join_audit, z_sha256=sha_file(self.z_path),
                    set_sha256=sha_file(self.set_path))

def split_genes(args, data):
    if args.smoke:
        order = np.random.default_rng(0).permutation(data.ng)
        test = np.sort(order[:data.ng // 5])
        train = np.sort(order[data.ng // 5:])
        order = np.random.default_rng(1).permutation(train)
        va = np.sort(order[:len(order) // 2])
        tr = np.sort(order[len(order) // 2:])
    else:
        test = np.array([i for i, g in enumerate(data.gl) if data.gchr[g] in FOLDS[args.fold]])
        train = np.array([i for i, g in enumerate(data.gl) if data.gchr[g] not in FOLDS[args.fold]])
        va_chr = set().union(*(FOLDS[(args.fold + j) % 5] for j in range(1, args.inner_va_folds + 1)))
        va = np.array([i for i in train if data.gchr[data.gl[i]] in va_chr])
        tr = np.array([i for i in train if data.gchr[data.gl[i]] not in va_chr])
    require(all(len(x) > 0 for x in (train, test, tr, va)), "empty chromosome split")
    require(not set(train) & set(test) and not set(tr) & set(va), "split overlap")
    require(set(tr) | set(va) == set(train), "inner partition mismatch")
    return train, test, tr, va

class Design:
    def __init__(self, args, data, gi):
        self.args = args
        self.arm = args.arm
        self.data = data
        self.train_pairs = np.flatnonzero(np.isin(data.gid, gi))
        self.is_train = np.zeros(len(data.keys), bool)
        self.is_train[self.train_pairs] = True
        self.basis = {}
        self.meta = {}
        self.blocks = []
        self.cluster = None
        for name in data.phi_continuous:
            x = data.col(name)
            mask = np.isfinite(x)
            if name == "dist_tss":
                require((x[mask] >= 0).all(), "negative dist_tss")
                x[mask] = np.log1p(x[mask])
            self.continuous(name, x, mask)
            self.meta[name]["na_pct"] = float(100 * (1 - mask.mean()))
        binaries = {name: data.col(name) for name in S3}
        for name in S3:
            require(np.isin(binaries[name], [0, 1]).all(), "invalid binary annotation")
        for name in S2:
            x = data.col(name)
            mask = np.isfinite(x)
            if name == "gh_link_score":
                require((x[mask] >= 0).all(), "negative link score")
                x[mask] = np.log1p(x[mask])
            self.continuous(name, x, mask)
            self.meta[name]["applicable_pct"] = float(100 * mask.mean())
            dup = [b for b in S3 if np.array_equal(binaries[b][self.train_pairs].astype(bool),
                                                  mask[self.train_pairs])]
            if dup:
                self.meta[name]["indicator"] = "dup of " + dup[0] + " — skipped"
            else:
                self.binary(name + "_A", mask.astype(float), "ind")
                self.meta[name]["indicator"] = "separate applicability indicator"
        for name in S3:
            self.binary(name, binaries[name])
        self.B = np.hstack([b[2] for b in self.blocks])
        self.names = [f"{b[0]}:{i}" for b in self.blocks for i in range(b[2].shape[1])]
        self.P = np.zeros((self.B.shape[1], self.B.shape[1]))
        self.block_penalties = []
        offset = 0
        for name, kind, block, penalty in self.blocks:
            width = block.shape[1]
            self.P[offset:offset + width, offset:offset + width] = penalty
            self.block_penalties.append((name, kind, penalty))
            offset += width
        self.blocks = None
        self.M = self.moment(self.B, self.train_pairs)
        self.initial = np.zeros(self.B.shape[1])
        if self.arm == "cluster":
            self.make_cluster()
        if self.arm == "nn":
            torch.manual_seed(0)
            self.net = torch.nn.Sequential(torch.nn.Linear(self.B.shape[1], 64),
                                          torch.nn.GELU(), torch.nn.Linear(64, 32),
                                          torch.nn.GELU(), torch.nn.Linear(32, 1)).to(args.device)
            torch.nn.init.zeros_(self.net[-1].weight)
            torch.nn.init.zeros_(self.net[-1].bias)
            self.initial = self.pack()
            self.nn_input_names = self.names.copy()
            self.names = [f"nn.{name}:{i}" for name, p in self.net.named_parameters()
                          for i in range(p.numel())]
            self.P = None
            self.M = None
        self.nparam = len(self.initial)
        require(np.isfinite(self.B).all(), "nonfinite design")
        if self.P is not None:
            require(np.linalg.eigvalsh(self.P).min() > 0, "full penalty is not positive definite")
        self.signature = digest(self.export())

    def moment(self, matrix, rows):
        result = np.zeros((matrix.shape[1], matrix.shape[1]))
        for lo in range(0, len(rows), self.args.design_chunk):
            block = matrix[rows[lo:lo + self.args.design_chunk]].astype(np.float64)
            result += block.T @ block
        return result / len(rows)

    def continuous(self, name, x, mask):
        tr = mask & self.is_train
        require(tr.any(), "annotation column has no observed fitting rows")
        observed = x[tr]
        mu = float(observed.mean()) if len(observed) else 0.0
        sd = float(observed.std()) + NUM_EPS if len(observed) else 1.0
        xs = np.zeros(len(x))
        xs[mask] = (x[mask] - mu) / sd
        self.basis[name] = dict(mu=mu, sd=sd, kind="lin")
        meta = dict(n_obs_train=int(tr.sum()))
        self.meta[name] = meta
        knots = np.quantile(xs[tr], np.linspace(0, 1, NKNOT)) if tr.any() else np.array([])
        unique = []
        for value in knots:
            if not unique or value - unique[-1] > 1e-6 * (knots[-1] - knots[0] + ZERO_SD):
                unique.append(value)
        if self.arm != "spline" or len(unique) < 3:
            size = float(np.square(xs[self.train_pairs]).mean())
            penalty = size if size > ZERO_SD else 1.0
            self.blocks.append((name, "lin", xs[:, None].astype(np.float32), np.array([[penalty]])))
            meta["fallback"] = "linear" if self.arm == "spline" else None
            meta["zero_column_penalty"] = bool(size <= ZERO_SD)
            return
        knots = np.array(unique)
        spl = SplineTransformer(degree=3, knots=knots[:, None], include_bias=True,
                                extrapolation="constant").fit(xs[tr, None])
        raw = spl.transform(xs[:, None])
        mean_tr = raw[tr].mean(axis=0)
        Q = null_space(np.ones((1, raw.shape[1])))
        block = ((raw - mean_tr) @ Q) * mask[:, None]
        gram = block[tr].T @ block[tr] / tr.sum()
        bs = spl.bsplines_[0]
        grid = np.linspace(bs.t[3], bs.t[-4], 2001)
        deriv = bs.derivative(2)(grid)
        omega = Q.T @ (deriv.T @ deriv * (grid[1] - grid[0])) @ Q
        ev = np.linalg.eigvalsh(omega)
        rank = int((ev > NUM_EPS * ev.max()).sum())
        require(rank > 0, "spline curvature rank zero")
        q = float(np.trace(np.linalg.solve(gram, omega)) / rank)
        require(np.isfinite(q) and q > 0, "invalid curvature scaling")
        penalty = gram + omega / q
        require(np.linalg.eigvalsh(penalty).min() > 0, "spline P not positive definite")
        self.blocks.append((name, "spl", block.astype(np.float32), penalty))
        self.basis[name].update(kind="spl", knots=knots.tolist(), mean_tr=mean_tr.tolist(), Q=Q.tolist())
        meta.update(rank_omega=rank, q=q, n_basis=raw.shape[1], knots=knots.tolist())

    def binary(self, name, values, kind="bin"):
        p = float(values[self.train_pairs].mean())
        self.blocks.append((name, kind, (values - p)[:, None].astype(np.float32), np.array([[1.0]])))
        self.basis[name] = dict(kind=kind, p=p)
        if name in S3:
            self.meta[name] = dict(p=p)

    def make_cluster(self):
        from sklearn.mixture import GaussianMixture
        settings = dict(covariance_type="diag", reg_covar=1e-6, n_init=3,
                        max_iter=200, tol=1e-3, random_state=20260907)
        grid = [4, 8, 16]
        train = np.ascontiguousarray(self.B[self.train_pairs])
        fitted = []
        for k in grid:
            model = GaussianMixture(n_components=k, **settings).fit(train)
            require(model.converged_, "GMM did not converge")
            fitted.append((float(model.bic(train)), k, model))
        _, k, model = min(fitted, key=lambda v: (v[0], v[1]))
        order = np.lexsort(tuple([np.arange(k)] + [model.means_[:, j] for j in range(self.B.shape[1] - 1, -1, -1)]))
        hard = self.args.cluster_assignment == "hard"
        membership = np.empty(len(self.B), np.int16) if hard else np.empty((len(self.B), k), np.float32)
        for lo in range(0, len(self.B), self.args.design_chunk):
            hi = min(lo + self.args.design_chunk, len(self.B))
            prob = model.predict_proba(self.B[lo:hi])[:, order]
            membership[lo:hi] = prob.argmax(axis=1) if hard else prob
        tm = membership[self.train_pairs]
        mass = np.bincount(tm, minlength=k).astype(float) if hard else tm.sum(axis=0, dtype=np.float64)
        require((mass > 0).all(), "empty selected cluster")
        weights = mass / mass.sum()
        Q = null_space(weights[None, :])
        require(np.max(np.abs(weights @ Q)) < ZERO_SD, "cluster weighted centering")
        self.P = Q.T @ np.diag(weights) @ Q
        self.B = Q[membership].astype(np.float32) if hard else (membership @ Q).astype(np.float32)
        self.M = self.moment(self.B, self.train_pairs)
        self.initial = np.zeros(k - 1)
        self.cluster = dict(selected_k=k, assignment=self.args.cluster_assignment,
                            k_grid=grid, bic={str(kk): bb for bb, kk, _ in fitted},
                            settings=settings, input_names=self.names,
                            train_mass_fraction=weights.tolist(), Q=Q.tolist(),
                            gmm_weights=model.weights_[order].tolist(),
                            gmm_means=model.means_[order].tolist(),
                            gmm_covariances=model.covariances_[order].tolist())
        self.names = [f"cluster_basis:{j}" for j in range(k - 1)]

    def pack(self):
        return torch.cat([p.detach().reshape(-1) for p in self.net.parameters()]).cpu().numpy().astype(np.float64)

    def unpack(self, a):
        offset = 0
        with torch.no_grad():
            for p in self.net.parameters():
                n = p.numel()
                p.copy_(torch.as_tensor(a[offset:offset + n], dtype=p.dtype,
                                        device=self.args.device).reshape(p.shape))
                offset += n

    def f_of(self, a, dtype=None):
        dtype = np.float32 if dtype is None else dtype
        result = np.empty(len(self.B), dtype=dtype)
        if self.arm == "nn":
            self.unpack(a)
        for lo in range(0, len(self.B), self.args.design_chunk):
            hi = min(lo + self.args.design_chunk, len(self.B))
            if self.arm == "nn":
                td = next(self.net.parameters()).dtype
                with torch.no_grad():
                    x = torch.as_tensor(self.B[lo:hi], device=self.args.device, dtype=td)
                    result[lo:hi] = self.net(x).squeeze(1).cpu().numpy()
            else:
                result[lo:hi] = self.B[lo:hi].astype(dtype, copy=False) @ a.astype(dtype)
        return result

    def backward(self, a, gm, dtype):
        if self.arm == "nn":
            self.unpack(a)
            self.net.zero_grad(set_to_none=True)
            td = next(self.net.parameters()).dtype
            for lo in range(0, len(self.B), self.args.design_chunk):
                hi = min(lo + self.args.design_chunk, len(self.B))
                x = torch.as_tensor(self.B[lo:hi], device=self.args.device, dtype=td)
                f = self.net(x).squeeze(1)
                f.backward(torch.as_tensor(gm[lo:hi], device=self.args.device, dtype=td))
            return torch.cat([p.grad.reshape(-1) for p in self.net.parameters()]).cpu().numpy().astype(np.float64)
        gradient = np.zeros(self.nparam)
        for lo in range(0, len(self.B), self.args.design_chunk):
            hi = min(lo + self.args.design_chunk, len(self.B))
            gradient += self.B[lo:hi].astype(dtype, copy=False).T @ gm[lo:hi]
        return gradient

    def penalty(self, a):
        pa = a if self.P is None else self.P @ a
        return float(a @ pa), 2 * pa

    def lambda_strong(self, gradient):
        if self.arm != "nn":
            pg = np.linalg.solve(self.P, gradient)
            return float(np.sqrt(max(0.0, pg @ self.M @ pg)) / (2 * EPS_F))
        self.unpack(self.initial)
        moment = np.zeros((33, 33))
        for lo in range(0, len(self.train_pairs), self.args.design_chunk):
            rows = self.train_pairs[lo:lo + self.args.design_chunk]
            with torch.no_grad():
                h = self.net[:-1](torch.as_tensor(self.B[rows], device=self.args.device))
                h = torch.cat([h, torch.ones((len(rows), 1), device=h.device)], axis=1)
                moment += (h.double().T @ h.double()).cpu().numpy()
        moment /= len(self.train_pairs)
        g = gradient[-33:]
        return float(np.sqrt(max(0.0, g @ moment @ g)) / (2 * EPS_F))

    def export(self):
        return dict(version="v9" if self.args.v9_enabled else "v8", basis=self.basis,
                    meta=self.meta, coef_names=self.names,
                    cols_S1=self.data.phi_continuous, cols_S2=S2, cols_S3=S3, CL=CL, NKNOT=NKNOT,
                    arm=self.arm, cluster=self.cluster, fold=self.args.fold,
                    smoke=self.args.smoke, band_root=str(self.data.root),
                    n_train_pairs=int(len(self.train_pairs)),
                    nn_architecture=[self.B.shape[1], 64, 32, 1] if self.arm == "nn" else None,
                    nn_input_names=self.nn_input_names if self.arm == "nn" else None,
                    transfer_note="linear/spline use exported basis; nn/cluster need their arm-specific reconstruction")

def mix_terms(T, pi, tau2, capped=False, cap=5.0):
    a = tau2 / (2 * (1 + tau2))
    log_lr = a * np.square(T) - 0.5 * np.log1p(tau2)
    signal = np.log(pi) + log_lr
    log_m = np.logaddexp(np.log1p(-pi), signal)
    posterior = np.exp(signal - log_m)
    derivative = -2 * a * T * posterior
    if capped:
        derivative = derivative * (log_m < cap)
        log_m = np.minimum(log_m, cap)
    return -log_m, derivative, posterior

def mix_em(T):
    pi, variance = 0.05, 4.0
    for _ in range(200):
        _, _, posterior = mix_terms(T, pi, variance - 1)
        new_pi = float(np.clip(posterior.mean(), 1e-3, 0.5))
        new_var = float(max(1.01, (posterior * T ** 2).sum() / max(posterior.sum(), ZERO_SD)))
        converged = abs(new_pi - pi) < 1e-7 and abs(new_var - variance) < 1e-6
        pi, variance = new_pi, new_var
        if converged:
            break
    return pi, variance - 1

def mixed_prior(pi_hat, pi_bbj, kappa):
    if kappa == 0.0:
        return np.broadcast_to(np.asarray(pi_hat, dtype=np.float64), np.shape(pi_bbj)).copy()
    if kappa == 1.0:
        return np.asarray(pi_bbj, dtype=np.float64).copy()
    return (1.0 - kappa) * np.asarray(pi_hat, dtype=np.float64) + kappa * np.asarray(pi_bbj, dtype=np.float64)

def target_prediction(T):
    return np.abs(T)

def zrow(S, device, dtype):
    Z = torch.as_tensor(S, device=device, dtype=dtype).clone()
    Z.sub_(Z.mean(axis=1, keepdim=True))
    sd = torch.sqrt(torch.mean(Z * Z, axis=1, keepdim=True) + NUM_EPS ** 2)
    Z.div_(sd)
    return Z, sd

class Trainer:
    def __init__(self, args, data, design, residuals, run_hash, truth=None):
        self.args, self.data, self.design = args, data, design
        self.residuals = residuals
        self.run_hash = run_hash
        self.truth = truth
        self.dtype = torch.float32
        self.R = torch.as_tensor(residuals.R, device=args.device)
        self.frozen = None
        self.gi = None
        self.kappa = 0.0

    def blocks(self, gi):
        for lo in range(0, len(gi), self.args.gene_batch):
            yield gi[lo:lo + self.args.gene_batch]

    def gene_T(self, phi, gi, conditioned=True):
        result = np.empty((len(gi), len(self.residuals.names)))
        offset = 0
        for block in self.blocks(gi):
            S = self.data.burden(phi * self.data.pw, block, conditioned=conditioned)
            Z, _ = zrow(S, self.args.device, self.dtype)
            T = (Z.double() @ self.R.double().T).cpu().numpy()
            result[offset:offset + len(block)] = T
            offset += len(block)
        require(np.isfinite(result).all(), "nonfinite T")
        return result

    def freeze_head(self, gi, existing=None):
        T = self.gene_T(np.ones(len(self.data.keys), np.float32), gi)
        existing = existing or {}
        return {t: existing[t] if t in existing else mix_em(T[:, i])
                for i, t in enumerate(self.residuals.names)}

    def head(self, S, gi, want_grad):
        Z, sd = zrow(S, self.args.device, self.dtype)
        T = (Z.double() @ self.R.double().T).cpu().numpy()
        dT = np.zeros_like(T)
        loss = 0.0
        if not self.args.v9_enabled:
            scale = len(T) * len(self.residuals.names)
            for it, trait in enumerate(self.residuals.names):
                terms, derivative, _ = mix_terms(T[:, it], *self.frozen[trait],
                                                 capped=self.args.obj == "mixbf_capped", cap=self.args.obj_cap)
                loss += float(terms.sum() / scale)
                dT[:, it] = derivative / scale
        else:
            require(self.truth is not None, "v9 trainer lacks BBJ truth")
            targets, _, pi_bbj, tested = self.truth.arrays(self.data, gi, self.residuals.names)
            if self.args.arm_prior_only:
                scale = len(T) * len(self.residuals.names)
                for it, trait in enumerate(self.residuals.names):
                    available_prior = np.where(np.isfinite(pi_bbj), pi_bbj, self.frozen[trait][0])
                    pi = mixed_prior(self.frozen[trait][0], available_prior, self.kappa)
                    terms, derivative, _ = mix_terms(T[:, it], pi, self.frozen[trait][1],
                                                     capped=self.args.obj == "mixbf_capped",
                                                     cap=self.args.obj_cap)
                    loss += float(terms.sum() / scale)
                    dT[:, it] = derivative / scale
            else:
                pi_hat = float(np.mean([self.frozen[t][0] for t in self.residuals.names]))
                available_prior = np.where(np.isfinite(pi_bbj), pi_bbj, pi_hat)
                weight = mixed_prior(pi_hat, available_prior, self.kappa)[:, None]
                prediction = target_prediction(T)
                delta = np.where(tested, prediction - targets, 0.0)
                loss = float(np.sum(weight * np.square(delta)))
                dT = 2.0 * weight * delta * np.sign(T)
        if not want_grad:
            return loss, None
        gZ = (torch.as_tensor(dT, device=self.args.device) @ self.R.double()).to(self.dtype)
        cross = (gZ * Z).mean(axis=1, keepdim=True)
        gZ.sub_(gZ.mean(axis=1, keepdim=True)).sub_(Z * cross).div_(sd)
        return loss, gZ.cpu().numpy()

    def data_loss(self, a, want_grad=True):
        dtype = np.float64 if self.dtype == torch.float64 else np.float32
        f = self.design.f_of(a, dtype)
        phi = np.exp(np.clip(f, -CL, CL)).astype(dtype)
        gradient_phi = np.zeros(len(f), dtype) if want_grad else None
        loss = 0.0
        for block in self.blocks(self.gi):
            S = self.data.burden(phi * self.data.pw, block)
            value, gS = self.head(S, block, want_grad)
            weight = 1.0 if self.args.v9_enabled and not self.args.arm_prior_only else len(block) / len(self.gi)
            loss += value * weight
            if want_grad:
                self.data.grad_phi(gS * weight, block, gradient_phi)
        if not want_grad:
            return loss, None
        gm = gradient_phi * self.data.pw * phi * (np.abs(f) < CL)
        return loss, self.design.backward(a, gm, dtype)

    def objective(self, a, lam):
        loss, gradient = self.data_loss(a)
        pen, gp = self.design.penalty(a)
        require(np.isfinite(loss) and np.isfinite(gradient).all(), "nonfinite objective/gradient")
        return loss + lam * pen, gradient + lam * gp

    def fit(self, lam, gi, frozen, maxiter, tag):
        self.gi, self.frozen = gi, frozen
        cfg = digest(dict(run=self.run_hash, design=self.design.signature, lam=lam,
                          gene_set=digest([self.data.gl[i] for i in gi]), frozen=frozen,
                          maxiter=maxiter, obj=self.args.obj, obj_cap=self.args.obj_cap,
                          kappa=self.kappa, arm_prior_only=self.args.arm_prior_only))
        root = Path(self.args.out) / tag
        ckpt = Path(str(root) + ".ckpt.npz")
        completed = Path(str(root) + ".fit.done")
        model_path = Path(str(root) + ".fit.a.npy")
        if completed.exists():
            with open(completed) as fh:
                saved = json.load(fh)
            require(saved["cfg"] == cfg, "completed fit configuration mismatch")
            require(saved["a_sha256"] == sha_file(model_path), "completed fit checksum mismatch")
            a = np.load(model_path, allow_pickle=False)
            require(a.shape == (self.design.nparam,) and np.isfinite(a).all(), "completed fit parameters")
            return a, saved["optimizer"]
        a = self.design.initial.copy()
        iteration = 0
        if ckpt.exists():
            with np.load(ckpt, allow_pickle=False) as z:
                require(str(z["cfg"]) == cfg, "checkpoint configuration mismatch")
                a = z["a"].astype(np.float64)
                iteration = int(z["iteration"])
            require(a.shape == self.design.initial.shape and np.isfinite(a).all(), "checkpoint parameters")
            log(f"resume {tag}, accepted iterations={iteration}; L-BFGS history restarts")
        require(iteration < maxiter, "exhausted incomplete checkpoint; explicit new budget required")
        count = [iteration]
        def checkpoint(x):
            tmp = Path(str(ckpt) + ".tmp")
            with open(tmp, "wb") as fh:
                np.savez(fh, a=x, iteration=count[0], cfg=cfg)
                fh.flush()
                os.fsync(fh.fileno())
            os.replace(tmp, ckpt)
        def callback(x):
            count[0] += 1
            if count[0] % 5 == 0:
                checkpoint(x)
                log(f"fit {tag}: accepted iterations={count[0]}")
        result = minimize(self.objective, a, args=(lam,), jac=True, method="L-BFGS-B",
                          callback=callback, options=dict(maxiter=maxiter - iteration,
                          maxfun=3 * (maxiter - iteration), ftol=1e-9, gtol=1e-10, maxcor=20))
        a = result.x.astype(np.float64)
        info = dict(nit=iteration + int(result.nit), nfev=int(result.nfev),
                    success=bool(result.success), status=int(result.status),
                    message=str(result.message), fun=float(result.fun),
                    resumed_from_iteration=iteration, resume_kind="parameter restart, not exact L-BFGS continuation")
        require(np.isfinite(a).all() and np.isfinite(result.fun), "invalid optimizer result")
        require(result.status in (0, 1), "optimizer failed beyond declared iteration/function budget")
        if self.args.require_convergence:
            require(result.success, "optimizer did not converge")
        atomic_npy(model_path, a)
        atomic_json(completed, dict(cfg=cfg, a_sha256=sha_file(model_path), optimizer=info))
        if ckpt.exists():
            ckpt.unlink()
        return a, info

    def gradcheck(self, gi, frozen):
        self.gi, self.frozen = gi, frozen
        self.dtype = torch.float64
        if self.design.arm == "nn":
            self.design.net.double()
        rng = np.random.default_rng(7)
        a = self.design.initial + rng.normal(0, 0.05, self.design.nparam)
        direction = rng.normal(size=self.design.nparam)
        direction /= np.linalg.norm(direction)
        step = 1e-3
        _, grad = self.objective(a, 1.0)
        plus, _ = self.objective(a + step * direction, 1.0)
        minus, _ = self.objective(a - step * direction, 1.0)
        fd, analytic = (plus - minus) / (2 * step), float(grad @ direction)
        relative = abs(fd - analytic) / max(abs(fd), abs(analytic), ZERO_SD)
        self.dtype = torch.float32
        if self.design.arm == "nn":
            self.design.net.float()
        require(relative < 0.05, "float64 directional gradient check")
        return dict(dtype="float64", finite_difference=fd, analytic=analytic, relative_error=relative)

def atomic_npy(path, value):
    path = owned_path(path)
    tmp = Path(str(path) + ".tmp")
    with open(tmp, "wb") as fh:
        np.save(fh, value, allow_pickle=False)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)

def metric_arrays(Tlearn, Tflat, names, frozen):
    flat = np.abs(Tflat).mean(axis=-2)
    require((flat > ZERO_SD).all(), "flat mean abs(T) zero; ratio undefined")
    ratio = np.abs(Tlearn).mean(axis=-2) / flat
    J = np.stack([mix_terms(Tlearn[..., i], *frozen[t])[0].mean(axis=-1)
                  for i, t in enumerate(names)], axis=-1)
    return ratio, J

def standardize_metric(observed, null, names, direction=1):
    observed = np.append(observed, np.mean(observed))
    values = np.column_stack([null, null.mean(axis=1)])
    mu, sd = values.mean(axis=0), values.std(axis=0, ddof=1)
    degenerate = sd <= ZERO_SD
    require(np.all(np.abs(observed[degenerate] - mu[degenerate]) <= NUM_EPS),
            "zero null SD with nonzero observed excess; cannot define z")
    scores = np.zeros_like(observed)
    null_z = np.zeros_like(values)
    good = ~degenerate
    scores[good] = direction * (observed[good] - mu[good]) / sd[good]
    null_z[:, good] = direction * (values[:, good] - mu[good]) / sd[good]
    result = {}
    for i, name in enumerate(names + ["mean"]):
        result[name] = dict(obs=float(observed[i]), null_mean=float(mu[i]), null_sd=float(sd[i]),
                            z=float(scores[i]), z_raw=float(direction * scores[i]),
                            z_direction="higher_better" if direction == 1 else "lower_J_is_positive_z",
                            degenerate_null=bool(degenerate[i]),
                            null_values=values[:, i].tolist())
    return result, null_z[:, -1]

def subset_ratio_stats(Tlearn, Tflat, mask, names):
    require(mask.shape == Tlearn.shape[1:] == Tflat.shape[1:], "U17 mask shape")
    require(mask.any(), "U17 scoring subset is empty")
    numerator = np.abs(Tlearn[:, mask]).mean(axis=1)
    denominator = np.abs(Tflat[:, mask]).mean(axis=1)
    require((denominator > ZERO_SD).all(), "U17 flat mean abs(T) zero")
    ratios = numerator / denominator
    observed, null = float(ratios[0]), ratios[1:]
    sd = float(null.std(ddof=1))
    if sd <= ZERO_SD:
        require(abs(observed - float(null.mean())) <= NUM_EPS,
                "U17 zero null SD with nonzero observed excess")
        z = 0.0
        null_z = np.zeros_like(null)
    else:
        z = float((observed - null.mean()) / sd)
        null_z = (null - null.mean()) / sd
    per_trait = {}
    for it, trait in enumerate(names):
        tm = mask[:, it]
        if not tm.any():
            per_trait[trait] = dict(n_genes=0, ratio=None)
            continue
        den = np.abs(Tflat[0, tm, it]).mean()
        per_trait[trait] = dict(n_genes=int(tm.sum()),
                                ratio=float(np.abs(Tlearn[0, tm, it]).mean() / den) if den > ZERO_SD else None)
    info = dict(ratio=observed, z=z, null_mean=float(null.mean()), null_sd=sd,
                q95=float(np.quantile(null, .95)),
                p=float((1 + np.sum(null >= observed)) / (len(null) + 1)),
                n_gene_trait_pairs=int(mask.sum()), n_unique_genes=int(mask.any(axis=1).sum()),
                per_trait=per_trait, null_values=null.tolist(),
                definition="mean abs(T_learned) / mean abs(T_flat) within the fixed subset")
    return info, null_z

def u17_metrics(Tlearn, Tflat, sig_mask, tested_mask, names):
    primary, primary_null_z = subset_ratio_stats(Tlearn, Tflat, sig_mask, names)
    secondary, _ = subset_ratio_stats(Tlearn, Tflat, tested_mask & ~sig_mask, names)
    passed = bool(primary["ratio"] > primary["q95"])
    primary["used_for_pass"] = True
    primary["pass"] = passed
    secondary["used_for_pass"] = False
    return dict(primary_bbj_designated=primary, secondary_bbj_nondesignated=secondary,
                pass_primary_only=passed,
                decision_rule="only primary ratio > its fixed-phi permutation q95"), primary_null_z

class Evaluation:
    def __init__(self, trainer, gi, frozen, count, seed, known, truth=None):
        self.trainer, self.gi, self.frozen = trainer, gi, frozen
        self.names = trainer.residuals.names
        self.known = known
        self.truth = truth
        self.perms = trainer.residuals.permutations(count, seed)
        h = hashlib.sha256()
        for perm in self.perms:
            h.update(perm.tobytes())
        self.plan = dict(B=count, seed=seed, permutation_sha256=h.hexdigest(),
                         conditioning="fixed phi and fixed train-flat mixture parameters",
                         exchangeability=trainer.residuals.perm_meta,
                         no_refitting=True, full_procedure_null=False)
        self.flat = self.scores(np.ones(len(trainer.data.keys), np.float32))

    def scores(self, phi):
        tr = self.trainer
        result = np.empty((1 + len(self.perms), len(self.gi), len(self.names)))
        offset = 0
        for block in tr.blocks(self.gi):
            S = tr.data.burden(phi * tr.data.pw, block)
            Z, _ = zrow(S, tr.args.device, tr.dtype)
            zd = Z.double()
            result[0, offset:offset + len(block)] = (zd @ tr.R.double().T).cpu().numpy()
            for lo in range(0, len(self.perms), tr.args.perm_batch):
                ps = self.perms[lo:lo + tr.args.perm_batch]
                matrix = np.concatenate([tr.residuals.R[:, p].T for p in ps], axis=1)
                rm = torch.as_tensor(matrix, dtype=torch.float64, device=tr.args.device)
                values = (zd @ rm).cpu().numpy().reshape(len(block), len(ps), len(self.names))
                result[1 + lo:1 + lo + len(ps), offset:offset + len(block)] = values.transpose(1, 0, 2)
            offset += len(block)
        require(np.isfinite(result).all(), "nonfinite observed/permuted T")
        return result

    def evaluate(self, phi=None):
        T = self.flat if phi is None else self.scores(phi)
        ratios, J = metric_arrays(T, self.flat, self.names, self.frozen)
        ratios_info, ratio_null_z = standardize_metric(ratios[0], ratios[1:], self.names)
        J_info, _ = standardize_metric(J[0], J[1:], self.names, direction=-1)
        flat_J = np.array([mix_terms(self.flat[0, :, i], *self.frozen[t])[0].mean()
                           for i, t in enumerate(self.names)])
        observed = dict(ratio={t: float(v) for t, v in zip(self.names, ratios[0])},
                        ratio_mean=float(ratios[0].mean()), J=float(J[0].mean()),
                        J_sum_legacy_scale=float(J[0].sum()),
                        J_flat=float(flat_J.mean()), J_per_trait={t: float(v) for t, v in zip(self.names, J[0])},
                        ratio_stats=ratios_info, J_stats=J_info,
                        retention=known_retention(self.known, self.trainer.data, self.gi,
                                                  self.names, T[0], self.flat[0]))
        null_ratio = ratios[1:]
        means = null_ratio.mean(axis=1)
        null = dict(self.plan, q95={t: float(np.quantile(null_ratio[:, i], .95)) for i, t in enumerate(self.names)},
                    q95_mean=float(np.quantile(means, .95)),
                    p={t: float((1 + np.sum(null_ratio[:, i] >= ratios[0, i])) / (len(means) + 1))
                       for i, t in enumerate(self.names)},
                    p_mean=float((1 + np.sum(means >= ratios[0].mean())) / (len(means) + 1)),
                    null_mean={t: float(null_ratio[:, i].mean()) for i, t in enumerate(self.names)},
                    null_mean_aggregate=float(means.mean()), null_sd_aggregate=float(means.std(ddof=1)))
        medflat = np.median(self.flat[0] ** 2, axis=0)
        inflation = np.divide(np.median(T[0] ** 2, axis=0), medflat,
                              out=np.full(len(self.names), np.nan), where=medflat > ZERO_SD)
        observed.update(holdout_learned={t: float(np.abs(T[0, :, i]).mean()) for i, t in enumerate(self.names)},
                        holdout_flat={t: float(np.abs(self.flat[0, :, i]).mean()) for i, t in enumerate(self.names)},
                        inflation_medT2={t: float(v) if np.isfinite(v) else None for t, v in zip(self.names, inflation)},
                        null_fixed_phi=null)
        if self.truth is not None:
            _, sig_mask, _, tested = self.truth.arrays(self.trainer.data, self.gi, self.names)
            observed["u17"], ratio_null_z = u17_metrics(T, self.flat, sig_mask,
                                                        tested, self.names)
        return observed, ratio_null_z

def read_known(path):
    if not path:
        return []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        require({"gene", "trait"} <= set(reader.fieldnames or []), "known-gene TSV schema")
        records = [(r["gene"].strip().split(".")[0], r["trait"].strip()) for r in reader]
    require(len(set(records)) == len(records), "duplicate known gene/trait key")
    require(all(t in LEGACY for _, t in records), "known genes limited to four legacy traits")
    return records

def known_retention(known, data, gi, names, learn, flat):
    loc = {data.gl[g].split(".")[0]: i for i, g in enumerate(gi)}
    all_genes = {g.split(".")[0] for g in data.gl}
    require(len(loc) == len(gi), "version-stripped gene IDs not unique")
    result = dict(enabled=bool(known), selection_use=False, per_trait={})
    for it, trait in enumerate(names):
        if trait not in LEGACY:
            continue
        requested = [g for g, t in known if t == trait]
        present = [g for g in requested if g in loc]
        usable = [g for g in present if abs(flat[loc[g], it]) > ZERO_SD]
        ratios = [float(abs(learn[loc[g], it]) / abs(flat[loc[g], it])) for g in usable]
        lost = [g for g, r in zip(usable, ratios) if r <= 0.5]
        result["per_trait"][trait] = dict(n_requested=len(requested), n_in_fold=len(present),
                                            n_evaluable=len(usable), n_zero_flat=len(present) - len(usable),
                                            retention=float(np.mean(ratios)) if ratios else None,
                                            n_drop_ge_50pct=len(lost), genes_drop_ge_50pct=lost,
                                            n_absent_from_cache=sum(g not in all_genes for g in requested),
                                            n_elsewhere_in_cache=sum(g in all_genes and g not in loc for g in requested),
                                            n_outside_fold_or_absent=len(requested) - len(present))
    return result

def select_one_se(candidates, null_z):
    keys = list(candidates)
    best = max(keys, key=lambda k: candidates[k]["ratio_stats"]["mean"]["z"])
    best_z = candidates[best]["ratio_stats"]["mean"]["z"]
    boundaries = {}
    eligible = []
    for key in keys:
        difference = null_z[best] - null_z[key]
        se = float(difference.std(ddof=1))
        gap = best_z - candidates[key]["ratio_stats"]["mean"]["z"]
        accepted = bool(gap <= se + ZERO_SD)
        boundaries[key] = dict(gap_to_best=float(gap), paired_null_se=se,
                               null_mean_mcse=se / np.sqrt(len(difference)),
                               minimum_eligible_z=float(best_z - se), eligible=accepted)
        if accepted:
            eligible.append(key)
    chosen = eligible[0]
    return chosen, dict(primary="ratio z of arithmetic mean of trait ratios", best=best,
                         best_z=float(best_z), chosen=chosen, strongest_first=keys,
                         paired_boundaries=boundaries, eligible=eligible,
                         flat_selected=chosen == "inf", flat_convention="z=0 for deterministic ratio 1",
                         se_formula="sd_b(z_best_null[b]-z_candidate_null[b]), ddof=1; no sqrt(B) divisor",
                         interpretation="fixed-phi paired-null 1-SE heuristic, not equivalence or CV sampling SE",
                         reason="strongest penalty within paired-null 1-SE of maximum ratio z")

def u17_primary_z(candidate):
    return candidate["u17"]["primary_bbj_designated"]["z"]

def select_one_se_u17(candidates, null_z):
    keys = list(candidates)
    best = max(keys, key=lambda key: u17_primary_z(candidates[key]))
    best_z = u17_primary_z(candidates[best])
    boundaries, eligible = {}, []
    for key in keys:
        difference = null_z[best] - null_z[key]
        se = float(difference.std(ddof=1))
        gap = best_z - u17_primary_z(candidates[key])
        accepted = bool(gap <= se + ZERO_SD)
        boundaries[key] = dict(gap_to_best=float(gap), paired_null_se=se,
                               null_mean_mcse=se / np.sqrt(len(difference)),
                               minimum_eligible_z=float(best_z - se), eligible=accepted)
        if accepted:
            eligible.append(key)
    chosen = eligible[0]
    return chosen, dict(primary="U17 BBJ-designated ratio permutation z", best=best,
                         best_z=float(best_z), chosen=chosen, strongest_first=keys,
                         paired_boundaries=boundaries, eligible=eligible,
                         flat_selected=chosen == "inf",
                         flat_convention="z=0 for deterministic ratio 1",
                         se_formula="sd_b(z_best_null[b]-z_candidate_null[b]), ddof=1; no sqrt(B) divisor",
                         interpretation="fixed-phi paired-null 1-SE heuristic",
                         reason="strongest penalty within paired-null 1-SE of maximum U17 z")

def select_kappa_one_se(candidates, null_z):
    kappas = sorted(candidates)
    best = max(kappas, key=lambda k: (u17_primary_z(candidates[k]), -k))
    best_z = u17_primary_z(candidates[best])
    eligible, rows = [], {}
    for kappa in kappas:
        difference = null_z[best] - null_z[kappa]
        se = float(difference.std(ddof=1))
        z = u17_primary_z(candidates[kappa])
        gap = best_z - z
        accepted = bool(gap <= se + ZERO_SD)
        if accepted:
            eligible.append(kappa)
        rows[f"{kappa:g}"] = dict(kappa=kappa, z=float(z), paired_null_se=se,
                                   gap_to_best=float(gap), eligible=accepted,
                                   selected=False)
    chosen = min(eligible)
    rows[f"{chosen:g}"]["selected"] = True
    return chosen, dict(candidates=rows, chosen=chosen, best=best, best_z=float(best_z),
                         eligible=eligible, tie_break="smaller kappa",
                         rule="each kappa at its fixed lambda_strong; choose smallest within paired-null 1-SE of best U17 z",
                         grid_rearranged_automatically=False,
                         se_formula="sd_b(z_best_null[b]-z_kappa_null[b]), ddof=1")

def kappa_eb_guide(trainer, gi, frozen, grid):
    T = trainer.gene_T(np.ones(len(trainer.data.keys), np.float32), gi)
    _, _, pi_bbj, _ = trainer.truth.arrays(trainer.data, gi, trainer.residuals.names)

    def value(kappa, derivatives=False):
        ll = score = curvature = 0.0
        for it, trait in enumerate(trainer.residuals.names):
            pi_hat, tau2 = frozen[trait]
            available_prior = np.where(np.isfinite(pi_bbj), pi_bbj, pi_hat)
            pi = mixed_prior(pi_hat, available_prior, float(kappa))
            terms, _, posterior = mix_terms(T[:, it], pi, tau2)
            ll -= float(terms.sum())
            if derivatives:
                delta = available_prior - pi_hat
                dlog_dpi = posterior / pi - (1.0 - posterior) / (1.0 - pi)
                score += float(np.sum(delta * dlog_dpi))
                curvature -= float(np.sum(np.square(delta * dlog_dpi)))
        return (ll, score, curvature) if derivatives else ll

    left = value(0.0, True)[1]
    right = value(1.0, True)[1]
    if left <= 0:
        khat = 0.0
    elif right >= 0:
        khat = 1.0
    else:
        lo, hi = 0.0, 1.0
        for _ in range(np.finfo(np.float64).nmant + 2):
            mid = (lo + hi) / 2.0
            if value(mid, True)[1] > 0:
                lo = mid
            else:
                hi = mid
        khat = (lo + hi) / 2.0
    ll, score, curvature = value(khat, True)
    return dict(training_fits=0, flat_only=True,
                curve=[dict(kappa=float(k), log_marginal=float(value(k))) for k in grid],
                kappa_hat=float(khat), log_marginal_hat=ll, score_at_hat=score,
                curvature_at_hat=curvature,
                grid_rearranged_automatically=False,
                instruction="guide only; a human must preregister any future grid change")

def staar_pair_weights(data, name):
    values = data.col(name)
    unique_keys, first, inverse = np.unique(data.keys, return_index=True, return_inverse=True)
    unique_values = values[first]
    mapped = unique_values[inverse]
    consistent = (np.isnan(values) & np.isnan(mapped)) | (values == mapped)
    require(consistent.all(), f"STAAR {name}: annotation differs across duplicate variant-gene rows")
    M = len(unique_keys)
    weights = np.full(M, 0.5, dtype=np.float64)
    observed = np.isfinite(unique_values)
    require(observed.any(), f"STAAR {name}: all annotation values missing")
    weights[observed] = rankdata(unique_values[observed], method="average") / M
    return weights[inverse].astype(np.float32), dict(M=M, n_observed=int(observed.sum()),
                                                     n_missing=int((~observed).sum()),
                                                     tie_method="average", missing_weight=0.5,
                                                     denominator="all unique band variants")

def acat_equivalent_z(scores):
    p = 2.0 * norm.sf(np.abs(scores))
    p = np.clip(p, np.nextafter(0.0, 1.0), np.nextafter(1.0, 0.0))
    cauchy = np.tan((0.5 - p) * np.pi).mean(axis=0)
    p_acat = 0.5 - np.arctan(cauchy) / np.pi
    p_acat = np.clip(p_acat, np.nextafter(0.0, 1.0), np.nextafter(1.0, 0.0))
    result = norm.isf(np.maximum(p_acat / 2.0, np.nextafter(0.0, 1.0)))
    require(np.isfinite(result).all(), "nonfinite STAAR ACAT equivalent z")
    return result

def evaluate_staar(args, trainer, gi, frozen, known, truth):
    evaluation = Evaluation(trainer, gi, frozen, args.null_perm, args.test_seed,
                            known, truth=truth)
    scores, per_annotation, rank_meta = [], {}, {}
    _, sig_mask, _, tested = truth.arrays(trainer.data, gi, trainer.residuals.names)
    for name in STAAR_ANNOTATIONS:
        phi, rank_meta[name] = staar_pair_weights(trainer.data, name)
        Tk = evaluation.scores(phi)
        scores.append(Tk)
        per_annotation[name], _ = u17_metrics(Tk, evaluation.flat, sig_mask,
                                               tested, trainer.residuals.names)
    combined = acat_equivalent_z(np.stack(scores, axis=0))
    combined_metrics, _ = u17_metrics(combined, evaluation.flat, sig_mask,
                                      tested, trainer.residuals.names)
    return dict(arm="staar", training_fits=0, maf_weight="Beta(1,1); no MAF weighting",
                annotations=STAAR_ANNOTATIONS, rank=rank_meta,
                gene_p="two-sided standard-normal p from the frozen v8 T score",
                combination="equal-weight ACAT (Cauchy)", permutations=evaluation.plan,
                per_annotation=per_annotation, u17=combined_metrics,
                flat_comparator="same genes, residuals, folds and permutations")

def verdicts(args, evaluation, clamp_fraction, correlations):
    names = list(evaluation["ratio"])
    ratio = dict(evaluation["ratio"], mean=evaluation["ratio_mean"])
    q95 = dict(evaluation["null_fixed_phi"]["q95"], mean=evaluation["null_fixed_phi"]["q95_mean"])
    infl = evaluation["inflation_medT2"].copy()
    infl["mean"] = float(np.mean(list(infl.values()))) if all(v is not None for v in infl.values()) else None
    first, second, third = {}, {}, {}
    for t in names + ["mean"]:
        first[t] = dict(ratio=ratio[t], q95=q95[t], pass_=bool(ratio[t] > q95[t]))
        bound = args.inflation_range
        second[t] = dict(value=infl[t], allowed=bound,
                          pass_=(bound[0] <= infl[t] <= bound[1]) if bound and infl[t] is not None else None)
        maximum = max(abs(v) for v in correlations.values())
        configured = args.clamp_max is not None and args.corr_max is not None
        third[t] = dict(clamp_frac=clamp_fraction, max_abs_corr=maximum,
                         clamp_max=args.clamp_max, corr_max=args.corr_max,
                         pass_=(clamp_fraction <= args.clamp_max and maximum <= args.corr_max) if configured else None,
                         shared_across_traits=True)
    return dict(i_ratio_vs_fixed_phi_q95=first, ii_inflation=second,
                iii_clamp_and_confounding=third,
                thresholds_status="configured" if args.inflation_range and args.clamp_max is not None and args.corr_max is not None else "pending requester thresholds; null means unjudged",
                interpretation="provisional diagnostics, not association p-values; trait-wise results unadjusted",
                smoke_only=bool(args.smoke), known_genes_are_not_a_gate=True)

def resolve_v9_inputs(args):
    explicit = [args.external_dir, args.bbj_gene_z, args.bbj_gene_set,
                args.cadd_side, args.cond_covariates, args.band]
    args.v9_enabled = bool(any(x is not None for x in explicit) or
                           args.arm == "staar" or args.arm_prior_only or
                           args.kappa_eb_guide)
    args.cond_cov_path = None
    args.cond_cov_manual_override = False
    if not args.v9_enabled:
        return
    require(args.band in ("ge05", "ge01"),
            "v9 requires --band ge05 (training) or ge01 (transfer)")
    directory = Path(args.external_dir) if args.external_dir else None
    args.bbj_gene_z = args.bbj_gene_z or (str(directory / "bbj_gene_z.tsv") if directory else None)
    args.bbj_gene_set = args.bbj_gene_set or (str(directory / "bbj_gene_set.tsv") if directory else None)
    args.cadd_side = args.cadd_side or (str(directory / "cadd_side.tsv") if directory else None)
    automatic = str(directory / f"cond_cov_{args.band}.npz") if directory else None
    if args.cond_covariates:
        args.cond_cov_path = args.cond_covariates
        args.cond_cov_manual_override = True
        log(f"WARNING manual conditioning override for --band {args.band}: automatic file selection bypassed")
    else:
        args.cond_cov_path = automatic
    require(all((args.bbj_gene_z, args.bbj_gene_set, args.cadd_side, args.cond_cov_path)),
            "v9 needs --external-dir or all explicit BBJ/CADD/conditioning paths")
    for label, path in (("bbj_gene_z", args.bbj_gene_z),
                        ("bbj_gene_set", args.bbj_gene_set),
                        ("cadd_side", args.cadd_side),
                        ("conditioning", args.cond_cov_path)):
        require(Path(path).is_file(), f"{label} input absent")

def validate_kappa_grid(grid):
    require(len(grid) >= 2, "kappa grid needs at least endpoints 0 and 1")
    require(all(math.isfinite(k) and 0 <= k <= 1 for k in grid),
            "kappa values must be finite in [0,1]")
    require(len(set(grid)) == len(grid), "duplicate kappa grid value")
    require(all(b > a for a, b in zip(grid, grid[1:])),
            "kappa grid must be strictly increasing")
    require(0.0 in grid and 1.0 in grid,
            "kappa grid must include both mandatory endpoints 0 and 1")

def conditional_diagnostics(args, data, residuals, path):
    require(data.cond is not None, "conditional diagnostics require conditioning data")
    names = [gene_key(x) for x in args.cond_check_genes]
    require(len(names) == len(set(names)), "duplicate conditional diagnostic gene")
    missing = [name for name in names if name not in data.gix_key]
    require(not missing, "conditional diagnostic genes absent from band: " + ",".join(missing))
    R = torch.as_tensor(residuals.R, device=args.device)
    known = {}
    for name in names:
        gi = np.asarray([data.gix_key[name]], dtype=np.int64)
        raw = data.burden(data.pw, gi, conditioned=False)
        cond = data.burden(data.pw, gi, conditioned=True)
        z_raw, _ = zrow(raw, args.device, torch.float32)
        z_cond, _ = zrow(cond, args.device, torch.float32)
        traw = (z_raw.double() @ R.double().T).cpu().numpy()[0]
        tcond = (z_cond.double() @ R.double().T).cpu().numpy()[0]
        known[name] = dict(raw_T=traw.tolist(), conditional_T=tcond.tolist(),
                           traits=residuals.names)
        log("COND_T " + json.dumps(dict(gene=name, raw_T=traw.tolist(),
                                         conditional_T=tcond.tolist()), allow_nan=False))
    lost = []
    for lo in range(0, data.ng, args.gene_batch):
        gi = np.arange(lo, min(lo + args.gene_batch, data.ng), dtype=np.int64)
        cond = data.burden(data.pw, gi, conditioned=True).astype(np.float64)
        cond_sd = cond.std(axis=1, ddof=0)
        dead = np.flatnonzero(cond_sd <= ZERO_SD)
        if len(dead):
            raw = data.burden(data.pw, gi[dead], conditioned=False).astype(np.float64)
            raw_sd = raw.std(axis=1, ddof=0)
            for local, rsd in zip(dead, raw_sd):
                lost.append(dict(gene_id=data.gl[int(gi[local])], raw_sd=float(rsd),
                                 conditional_sd=float(cond_sd[local]),
                                 reason="conditional burden has numerical zero variance"))
    atomic_tsv(path, ["gene_id", "raw_sd", "conditional_sd", "reason"], lost)
    return dict(known_gene_T=known, configured_genes=names, all_configured_present=True,
                lost_gene_count=len(lost), lost_definition=f"conditional burden population SD <= {ZERO_SD:g}",
                lost_genes_are_diagnostic_not_discoveries=True,
                inverse_calculations=data.cond["inverse_calculations"],
                covariates=data.cond["n_covariates"], id_join=data.cond["id_join"],
                snp_list_sha256=data.cond["snp_list_sha256"])

def parser():
    ap = argparse.ArgumentParser(description=__doc__)
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--fold", type=int, choices=range(5))
    mode.add_argument("--smoke", choices=[str(i) for i in range(1, 23)])
    mode.add_argument("--make-smoke-fixture", metavar="DIR")
    ap.add_argument("--arm", default="spline", choices=["linear", "spline", "nn", "cluster", "staar"])
    ap.add_argument("--root", default=str(REF))
    ap.add_argument("--out", default=str(REF / "annot/l1"))
    ap.add_argument("--traits", default=str(PREREG))
    ap.add_argument("--resid", choices=["eta", "cov"], default="cov")
    ap.add_argument("--resid-dir", default=str(OWNER / "work/ref/annot/resid"))
    ap.add_argument("--offset-dir", default=str(REF / "annot/offset"))
    ap.add_argument("--holdout-trait", action="append", default=[])
    ap.add_argument("--known-genes")
    ap.add_argument("--external-dir",
                    help="directory containing exact E1/E3/E4 filenames from the v9 brief")
    ap.add_argument("--bbj-gene-z", help="manual bbj_gene_z.tsv path")
    ap.add_argument("--bbj-gene-set", help="manual bbj_gene_set.tsv path")
    ap.add_argument("--cadd-side", help="manual cadd_side.tsv path")
    ap.add_argument("--band", choices=["ge05", "ge01"],
                    help="ge05 selects training conditioning; ge01 selects transfer conditioning")
    ap.add_argument("--cond-covariates",
                    help="manual conditioning NPZ override; emits a warning because --band normally selects it")
    ap.add_argument("--target-clip-q", type=float, default=0.99)
    ap.add_argument("--bbj-prior-sigmoid-scale", type=float, default=1.0,
                    help="scale applied to standardized abs(z_cond) inside sigmoid; 1 is the stated sigmoid")
    ap.add_argument("--kappa-grid", type=float, nargs="+", default=DEFAULT_KAPPA_GRID)
    ap.add_argument("--kappa-eb-guide", action="store_true",
                    help="write flat-only marginal-likelihood guide and exit without phi fitting")
    ap.add_argument("--arm-prior-only", action="store_true",
                    help="S1 lower-bound control: retain v8 phenotype loss and replace only its gene prior")
    ap.add_argument("--min-bbj-genes", type=int, default=100)
    ap.add_argument("--cond-check-genes", nargs="+", default=DEFAULT_COND_CHECK_GENES)
    ap.add_argument("--synthetic-smoke", action="store_true",
                    help="fixture-only relaxation of the fixed cohort sample-count gate")
    ap.add_argument("--perm-strata", help="optional IID,stratum TSV, e.g. cohort/ancestry blocks")
    ap.add_argument("--perm-r", type=int, default=0)
    ap.add_argument("--cluster-assignment", choices=["hard", "soft"], default="hard")
    ap.add_argument("--export-basis-only", action="store_true")
    ap.add_argument("--gradcheck-only", action="store_true")
    ap.add_argument("--require-convergence", action="store_true")
    ap.add_argument("--obj", choices=["mixbf", "mixbf_capped"], default="mixbf")
    ap.add_argument("--obj-cap", type=float, default=5.0)
    ap.add_argument("--inner-va-folds", type=int, choices=[1, 2], default=2)
    ap.add_argument("--lam", type=float, help="fixed lambda sensitivity only; inf selects flat")
    ap.add_argument("--learning-curve", type=float, nargs="+", default=None,
                    help="[병목 진단 2026-09-08] inner_tr 유전자 비율 목록; 각 비율×격자 λ 로 내부 적합 후 같은 inner_va·순열로 z 기록하고 종료(최종 적합 없음). 설계·매듭·head 는 전체 inner_tr 고정.")
    ap.add_argument("--lam-divisors", type=float, nargs="+", default=[10.0, 100.0, 1000.0],
                    help="lambda_strong divisors, strongest first; default D9 grid. Fix before viewing outcomes.")
    ap.add_argument("--maxiter", type=int, default=60)
    ap.add_argument("--inner-maxiter", type=int, default=30)
    ap.add_argument("--inner-null-perm", type=int, default=50)
    ap.add_argument("--null-perm", type=int, default=100)
    ap.add_argument("--inner-seed", type=int, default=20260908)
    ap.add_argument("--test-seed", type=int, default=123)
    ap.add_argument("--gene-batch", type=int, default=256)
    ap.add_argument("--perm-batch", type=int, default=4)
    ap.add_argument("--design-chunk", type=int, default=16384)
    ap.add_argument("--min-n", type=int, default=1000)
    ap.add_argument("--inflation-range", nargs=2, type=float)
    ap.add_argument("--clamp-max", type=float)
    ap.add_argument("--corr-max", type=float)
    safety_args(ap)
    return ap

def export_basis(path, design, code_sha):
    import pickle
    value = design.export()
    value["code_sha"] = code_sha
    tmp = Path(str(path) + ".tmp")
    with open(tmp, "wb") as fh:
        pickle.dump(value, fh, protocol=pickle.HIGHEST_PROTOCOL)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)

def make_smoke_fixture(target):
    global np
    import numpy as np
    import scipy.sparse as sparse

    target = owned_path(target)
    try:
        target.relative_to(Path.cwd().resolve())
    except ValueError:
        require(False, "synthetic fixture must be inside the current work directory")
    require(not target.exists(), "synthetic fixture target already exists")
    root = target / "root"
    external = target / "external"
    resid = target / "resid"
    output = target / "out"
    for directory in (root / "annot/cache", root / "annot/ds", root / "groupfiles_bwg",
                      external, resid, output):
        directory.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(20260908)
    n_samples, n_genes, variants_per_gene, n_traits = 160, 600, 3, 25
    lipid = ["LDLR", "HMGCR", "PCSK9", "APOB", "LPA"]
    genes = lipid + [f"ENSGSYN{i:09d}" for i in range(n_genes - len(lipid))]
    keys = np.asarray([f"synthetic_variant_{i:06d}" for i in range(n_genes * variants_per_gene)])
    pair_genes = np.repeat(np.asarray(genes), variants_per_gene)
    pairs = len(keys)
    required = S1 + S2 + S3 + DIAG + ["n_genes"]
    extras = [f"unused_{i:02d}" for i in range(33 - len(required))]
    cols = required + extras
    require(len(cols) == 33 and len(set(cols)) == 33, "internal smoke cache column construction")
    Xann = rng.normal(size=(pairs, 33)).astype(np.float32)
    for name in ("dist_tss", "gh_link_score", "maf", "r2", "avg_cs"):
        Xann[:, cols.index(name)] = np.abs(Xann[:, cols.index(name)])
    for name in S3:
        Xann[:, cols.index(name)] = rng.integers(0, 2, pairs)
    t1_na = rng.random(pairs) < 0.05
    Xann[:, cols.index("t1_na")] = t1_na
    Xann[:, cols.index("n_genes")] = 1.0
    for j, name in enumerate(S2):
        missing = (np.arange(pairs) + j) % (11 + j) == 0
        Xann[missing, cols.index(name)] = np.nan
    lipid_pairs = np.isin(pair_genes, lipid)
    for name in ("cons", "linsight", "gpn_msa"):
        Xann[lipid_pairs, cols.index(name)] += 3.0
    atomic_npz(root / "annot/cache/fm_all.npz", X=Xann,
               cols=np.asarray(cols), key37=keys, gene=pair_genes,
               chr=np.full(pairs, 22, dtype=np.int16))
    atomic_json(root / "groupfiles_bwg/chr22.summary.json", {"pairs": pairs})

    samples = np.asarray([f"SYNTH_SAMPLE_{i:04d}" for i in range(n_samples)])
    Xcov = rng.binomial(2, 0.25, size=(n_samples, 5)).astype(np.float64)
    dosage = rng.binomial(2, 0.06, size=(pairs, n_samples)).astype(np.float32)
    linked = rng.random(pairs) < 0.35
    dosage[linked] = np.minimum(2.0, dosage[linked] + Xcov[:, 0])
    dosage[lipid_pairs] = np.minimum(2.0, dosage[lipid_pairs] + Xcov[:, 1])
    csr = sparse.csr_matrix(dosage, dtype=np.float32)
    atomic_npz(root / "annot/ds/chr22.ds.npz", samples=samples,
               indptr=csr.indptr.astype(np.int32), indices=csr.indices.astype(np.int32),
               data=csr.data.astype(np.float32), keys=keys)

    shuffled = rng.permutation(n_samples)
    cond_arrays = dict(ids=samples[shuffled], X=Xcov[shuffled],
                       snp_list=np.asarray([f"synthetic_lead_{i}" for i in range(5)]))
    atomic_npz(external / "cond_cov_ge05.npz", **cond_arrays)
    atomic_npz(external / "cond_cov_ge01.npz", **cond_arrays)
    cadd = np.abs(rng.normal(15.0, 5.0, pairs))
    cadd[lipid_pairs] += 15.0
    atomic_tsv(external / "cadd_side.tsv", ["variant_key", "cadd_phred"],
               [dict(variant_key=key, cadd_phred=f"{value:.12g}")
                for key, value, missing in zip(keys, cadd, t1_na) if not missing])

    traits = [f"trait{i:02d}" for i in range(1, n_traits + 1)]
    traits_path = target / "traits.tsv"
    atomic_tsv(traits_path, ["trait", "type", "transform", "file", "col"],
               [dict(trait=t, type="quant", transform="none", file="synthetic", col=t)
                for t in traits])
    z_rows, set_rows, set_50_rows = [], [], []
    designated = set(genes[::2])
    designated_50 = set(genes[:50])
    for gene in genes:
        for trait in traits:
            z = float(rng.normal() + (6.0 if gene in lipid else 0.0))
            z_rows.append(dict(gene_id=gene, trait=trait, z_cond=f"{z:.12g}",
                               n_snps=variants_per_gene))
            set_rows.append(dict(gene_id=gene, trait=trait,
                                 is_bbj_sig=int(gene in designated)))
            set_50_rows.append(dict(gene_id=gene, trait=trait,
                                    is_bbj_sig=int(gene in designated_50)))
    atomic_tsv(external / "bbj_gene_z.tsv", ["gene_id", "trait", "z_cond", "n_snps"], z_rows)
    atomic_tsv(external / "bbj_gene_set.tsv", ["gene_id", "trait", "is_bbj_sig"], set_rows)
    atomic_tsv(external / "bbj_gene_set_50.tsv", ["gene_id", "trait", "is_bbj_sig"], set_50_rows)

    residual_hashes = {}
    for it, trait in enumerate(traits):
        residual = (rng.normal(size=n_samples) + 0.25 * Xcov[:, it % Xcov.shape[1]] +
                    0.10 * dosage[(it * variants_per_gene) % pairs]).astype(np.float64)
        path = resid / f"{trait}.resid.tsv"
        atomic_tsv(path, ["IID", "resid"],
                   [dict(IID=iid, resid=f"{value:.12g}") for iid, value in zip(samples, residual)])
        residual_hashes[trait] = {"residual_sha256": sha_file(path)}
    summary_path = resid / "build_resid_v8.summary.json"
    atomic_json(summary_path, dict(traits_sha256=sha_file(traits_path),
                                   sample_order_sha256=digest(samples.tolist()),
                                   traits=residual_hashes))
    atomic_json(resid / "build_resid_v8.done", {"summary_sha256": sha_file(summary_path)})
    atomic_json(target / "fixture_manifest.json",
                dict(synthetic=True, chromosome="22", samples=n_samples, genes=n_genes,
                     variants=pairs, traits=n_traits, conditioning_columns=5,
                     bbj_designated_genes=len(designated), boosted_lipid_genes=len(lipid),
                     contains_real_participant_or_variant_values=False,
                     root=str(root), external_dir=str(external), resid_dir=str(resid),
                     traits_path=str(traits_path), out=str(output)))
    print("SMOKE_FIXTURE_DONE " + str(target), flush=True)

def smoke_contract_checks():
    import importlib.util

    T = np.asarray([-2.0, -0.25, 0.5, 1.75], dtype=np.float64)
    pi_hat, tau2 = 0.17, 3.0
    pi_bbj = np.asarray([0.1, 0.3, 0.7, 0.9], dtype=np.float64)
    old_terms, old_derivative, _ = mix_terms(T, pi_hat, tau2)
    new_terms, new_derivative, _ = mix_terms(T, mixed_prior(pi_hat, pi_bbj, 0.0), tau2)
    require(np.array_equal(old_terms, new_terms) and np.array_equal(old_derivative, new_derivative),
            "smoke regression: kappa=0 prior-only differs numerically from v8 phenotype loss")
    endpoint = mixed_prior(pi_hat, pi_bbj, 1.0)
    require(np.array_equal(endpoint, pi_bbj),
            "smoke regression: kappa=1 prior differs from BBJ prior")
    spec = importlib.util.spec_from_file_location("l1_train_v8_regression", Path(__file__).with_name("l1_train_v8.py"))
    v8 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(v8)
    v8.np, v8.torch = np, torch
    names = ["trait01", "trait02"]
    R = np.asarray([[0.2, -0.1, 0.3, -0.4, 0.15, -0.05],
                    [-0.3, 0.25, 0.1, -0.2, 0.05, 0.1]], dtype=np.float32)
    S = np.asarray([[0.0, 1.0, 0.5, 2.0, -0.5, 0.25],
                    [1.5, 0.2, -0.1, 0.8, 0.4, -0.7],
                    [0.3, -0.8, 1.1, 0.6, -0.2, 0.9]], dtype=np.float32)
    frozen = {"trait01": (0.17, 3.0), "trait02": (0.23, 2.0)}
    ns = argparse.Namespace
    old = object.__new__(v8.Trainer)
    old.args = ns(device="cpu", obj="mixbf", obj_cap=5.0)
    old.residuals = ns(R=R, names=names)
    old.R = torch.as_tensor(R)
    old.dtype, old.frozen = torch.float32, frozen
    old_loss, old_grad = old.head(S, True)

    class ContractTruth:
        @staticmethod
        def arrays(data, gi, trait_names):
            return (np.zeros((len(gi), len(trait_names))),
                    np.zeros((len(gi), len(trait_names)), dtype=bool),
                    np.asarray([0.1, 0.4, 0.8], dtype=np.float64)[:len(gi)],
                    np.ones((len(gi), len(trait_names)), dtype=bool))

    new = object.__new__(Trainer)
    new.args = ns(device="cpu", obj="mixbf", obj_cap=5.0,
                  v9_enabled=True, arm_prior_only=True)
    new.data, new.truth = None, ContractTruth()
    new.residuals = ns(R=R, names=names)
    new.R = torch.as_tensor(R)
    new.dtype, new.frozen, new.kappa = torch.float32, frozen, 0.0
    new_loss, new_grad = new.head(S, np.arange(len(S)), True)
    loss_diff = abs(old_loss - new_loss)
    grad_diff = float(np.max(np.abs(old_grad - new_grad)))
    require(loss_diff == 0.0 and grad_diff == 0.0,
            "smoke regression: kappa=0 prior-only head differs from imported v8")

    class ContractData:
        def __init__(self):
            self.keys = np.asarray([f"v{i}" for i in range(6)])
            self.pw = np.ones(6, dtype=np.float32)
            self.pairs_of = [np.asarray([0, 1]), np.asarray([2, 3]), np.asarray([4, 5])]
            self.D = np.asarray([[0.0, 1.0, 0.5, 0.0, 2.0, 0.2],
                                 [1.0, 0.0, 0.2, 1.5, 0.0, 0.4],
                                 [0.1, 0.5, 1.0, 0.0, 0.3, 1.2],
                                 [0.0, 1.1, 0.0, 0.6, 0.8, 0.0],
                                 [0.7, 0.0, 1.3, 0.2, 0.0, 0.5],
                                 [0.2, 0.9, 0.0, 1.0, 0.4, 0.1]], dtype=np.float32)

        def burden(self, weights, gi, conditioned=True):
            return np.asarray([weights[self.pairs_of[g]] @ self.D[self.pairs_of[g]]
                               for g in gi], dtype=weights.dtype)

        def grad_phi(self, gS, gi, result, conditioned=True):
            for row, g in enumerate(gi):
                pp = self.pairs_of[g]
                result[pp] = self.D[pp] @ gS[row]

    class ContractDesign:
        def __init__(self):
            self.B = np.asarray([[0.2, -0.1], [0.1, 0.3], [-0.2, 0.4],
                                 [0.5, 0.2], [-0.3, -0.1], [0.4, -0.2]], dtype=np.float32)

        def f_of(self, a, dtype=None):
            return self.B.astype(dtype) @ np.asarray(a, dtype=dtype)

        def backward(self, a, gm, dtype):
            return self.B.astype(dtype).T @ gm

    contract_data, contract_design = ContractData(), ContractDesign()
    common_args = dict(device="cpu", obj="mixbf", obj_cap=5.0, gene_batch=2,
                       arm_prior_only=True)
    old.args = ns(**common_args, v9_enabled=False)
    old.data, old.design, old.gi = contract_data, contract_design, np.arange(3)
    new.args = ns(**common_args, v9_enabled=True)
    new.data, new.design, new.gi = contract_data, contract_design, np.arange(3)
    a_contract = np.asarray([0.15, -0.2], dtype=np.float64)
    old_total, old_total_grad = old.data_loss(a_contract)
    new_total, new_total_grad = new.data_loss(a_contract)
    total_loss_diff = abs(old_total - new_total)
    total_grad_diff = float(np.max(np.abs(old_total_grad - new_total_grad)))
    require(total_loss_diff == 0.0 and total_grad_diff == 0.0,
            "smoke regression: kappa=0 prior-only data_loss differs from imported v8")
    missing_endpoint_failed = False
    try:
        validate_kappa_grid([0.0, 0.25, 0.5, 0.75])
    except RuntimeError as error:
        missing_endpoint_failed = str(error).startswith("GATE FAIL:")
    require(missing_endpoint_failed, "smoke regression: missing-kappa endpoint gate did not fail")
    log("SMOKE_CONTRACT imported_v8_data_loss_max_abs_diff=0 "
        "imported_v8_data_gradient_max_abs_diff=0 kappa1_bbj_max_abs_diff=0 "
        "missing_kappa_1_gate=PASS")

def verify_input_fingerprints(data):
    for name, before in data.fingerprints.items():
        st = Path(name).stat()
        require(before == dict(size=st.st_size, mtime_ns=st.st_mtime_ns),
                "input changed during run")

def main():
    args = parser().parse_args()
    if args.make_smoke_fixture:
        make_smoke_fixture(args.make_smoke_fixture)
        return
    resolve_v9_inputs(args)
    require(not args.gradcheck_only or args.smoke, "gradcheck-only requires smoke")
    require(not (args.gradcheck_only and args.export_basis_only), "incompatible diagnostic modes")
    require(not args.synthetic_smoke or args.smoke, "--synthetic-smoke requires --smoke")
    require(args.inner_null_perm >= 2 and args.null_perm >= 2, "at least two null permutations required")
    require(args.inner_seed != args.test_seed, "inner/test permutation seeds must differ")
    require(args.lam is None or (args.lam > 0 and not math.isnan(args.lam)), "lambda must be positive or inf")
    _dv = list(args.lam_divisors)
    require(all(math.isfinite(v) and v > 0 for v in _dv), "divisors must be positive and finite")
    require(all(b > a for a, b in zip(_dv, _dv[1:])), "divisors must be strictly increasing (strongest first)")
    _keys = ["div" + f"{v:g}" for v in _dv]
    require(len(set(_keys)) == len(_keys), "divisor keys collide at %g precision: " + ",".join(_keys))
    args.custom_grid = _dv != [10.0, 100.0, 1000.0]
    require(args.obj_cap > 0 and math.isfinite(args.obj_cap), "cap must be positive finite")
    require(all(v > 0 for v in (args.maxiter, args.inner_maxiter, args.min_n,
                                args.gene_batch, args.perm_batch, args.design_chunk)), "positive budgets required")
    require(args.inflation_range is None or 0 < args.inflation_range[0] < args.inflation_range[1], "inflation bounds")
    require(args.clamp_max is None or 0 <= args.clamp_max <= 1, "clamp bound")
    require(args.corr_max is None or 0 <= args.corr_max <= 1, "correlation bound")
    if args.v9_enabled:
        validate_kappa_grid(list(args.kappa_grid))
        require(0 < args.target_clip_q <= 1, "target clip quantile must be in (0,1]")
        require(args.min_bbj_genes > 0, "minimum BBJ gene count must be positive")
        require(math.isfinite(args.bbj_prior_sigmoid_scale) and args.bbj_prior_sigmoid_scale > 0,
                "BBJ prior sigmoid scale must be positive finite")
        require(not (args.arm == "staar" and args.arm_prior_only),
                "STAAR and prior-only arms are mutually exclusive")
        require(not (args.arm == "staar" and (args.learning_curve or args.gradcheck_only or
                                               args.export_basis_only or args.lam is not None)),
                "STAAR has zero training and cannot use fitting/basis modes")
        require(not (args.kappa_eb_guide and args.arm == "staar"),
                "kappa EB guide and STAAR are separate zero-training modes")
        if args.smoke:
            require(args.smoke == "22", "v9 smoke is restricted to chromosome 22")
            require(args.threads == 4, "v9 smoke requires exactly 4 threads")
            require(args.memory_gb <= 8, "v9 smoke memory cap must be <=8 GiB")
        if args.band == "ge01":
            require(args.arm == "staar" or args.export_basis_only or args.gradcheck_only,
                    "transfer band ge01 cannot drive kappa/lambda/model selection")
    rows, schema = read_traits(args.traits)
    names = [r["trait"] for r in rows]
    require(set(args.holdout_trait) <= set(names), "unknown holdout trait")
    active_rows = [r for r in rows if r["trait"] not in args.holdout_trait]
    require(bool(active_rows), "all traits held out")
    code_sha = sha_file(__file__)
    lock, guard = resource_guard(args)
    load_libraries(args.threads)
    require(args.device != "cuda" or torch.cuda.is_available(), "requested CUDA unavailable")
    data = Data(args)
    train, test, inner_tr, inner_va = split_genes(args, data)
    truth = BBJTruth(args, names) if args.v9_enabled else None
    if truth is not None:
        truth.validate_against_data(data, [r["trait"] for r in active_rows], args.min_bbj_genes)
        if args.synthetic_smoke:
            smoke_contract_checks()
    config = dict(vars(args))
    config["lam"] = "inf" if args.lam == float("inf") else args.lam
    manifest = dict(config=config, code_sha256=code_sha, traits_sha256=sha_file(args.traits),
                    traits_schema=schema, traits=rows, input_stat_fingerprints=data.fingerprints,
                    known_genes_sha256=sha_file(args.known_genes) if args.known_genes else None,
                    perm_strata_sha256=sha_file(args.perm_strata) if args.perm_strata else None)
    if truth is not None:
        manifest.update(version="v9", bbj=truth.manifest(), cadd_join=data.cadd_join,
                        conditioning=dict(band=args.band,
                                          automatic_file=f"cond_cov_{args.band}.npz",
                                          selected_path=args.cond_cov_path,
                                          manual_override=args.cond_cov_manual_override,
                                          n_covariates=data.cond["n_covariates"] if data.cond else None,
                                          id_join=data.cond["id_join"] if data.cond else None,
                                          inverse_calculations=data.cond["inverse_calculations"] if data.cond else None))
    if args.export_basis_only:
        design = Design(args, data, train)
        version = "v9" if args.v9_enabled else "v8"
        tag = f"basis_{args.arm}_{version}_" + (f"smoke_chr{args.smoke}" if args.smoke else f"fold{args.fold}")
        path = Path(args.out) / (tag + ".pkl")
        export_basis(path, design, code_sha)
        atomic_json(Path(args.out) / (tag + ".done"), dict(basis_sha256=sha_file(path), manifest=manifest))
        print("L1_DONE basis_export", flush=True)
        lock.close()
        return
    residuals = Residuals(args, data, active_rows)
    manifest["active_residual_sha256"] = residuals.hashes
    if args.resid == "cov":
        manifest["residual_bundle_sha256"] = sha_file(Path(args.resid_dir) / "build_resid_v8.summary.json")
    else:
        manifest["residual_bundle_sha256"] = digest({t: sha_file(Path(args.offset_dir) / f"{t}.eta.tsv") for t in names})
    manifest["resid_kind"] = residuals.kind
    manifest["sample_order_sha256"] = digest(data.samples)
    run_hash = digest(manifest)
    version = "v9" if args.v9_enabled else "v8"
    objective_arm = "prior_only" if args.arm_prior_only else args.arm
    tag = f"{objective_arm}_{version}_" + (f"smoke_chr{args.smoke}" if args.smoke else f"fold{args.fold}")
    tag += "_" + run_hash[:16]
    final_marker = Path(args.out) / (tag + ".done")
    if final_marker.exists():
        with open(final_marker) as fh:
            done = json.load(fh)
        require(done["run_hash"] == run_hash and done["status"] == "complete", "run completion configuration")
        output_root = Path(args.out).resolve()
        for name, expected in done["artifacts"].items():
            artifact = (output_root / name).resolve()
            try:
                artifact.relative_to(output_root)
            except ValueError:
                require(False, "run completion artifact escapes output directory")
            require(artifact.is_file() and sha_file(artifact) == expected,
                    "run completion artifact checksum")
        log("RESULT verified existing completed run " + tag)
        print("L1_DONE", flush=True)
        lock.close()
        return
    manifest_path = Path(args.out) / (tag + ".manifest.json")
    atomic_json(manifest_path, manifest)
    run_dir = Path(args.out) / tag
    if args.v9_enabled:
        run_dir.mkdir(parents=True, exist_ok=True)
    known = read_known(args.known_genes)
    cond_diagnostics = None
    cond_lost_path = None
    if args.v9_enabled:
        cond_lost_path = run_dir / "cond_lost_genes.tsv"
        cond_diagnostics = conditional_diagnostics(args, data, residuals, cond_lost_path)

    if args.arm == "staar":
        trainer = Trainer(args, data, None, residuals, run_hash, truth=truth)
        frozen = trainer.freeze_head(train)
        staar_result = evaluate_staar(args, trainer, test, frozen, known, truth)
        evaluation_path = run_dir / "evaluation.json"
        atomic_json(evaluation_path, staar_result)
        result_path = Path(args.out) / (tag + ".json")
        atomic_json(result_path, dict(version="v9", code_sha256=code_sha, run_hash=run_hash,
                                      mode="smoke" if args.smoke else f"fold{args.fold}",
                                      arm="staar", evaluation=staar_result,
                                      conditional=cond_diagnostics, bbj=truth.manifest(),
                                      resource_guard=guard, sec=round(time.time() - T0)))
        artifacts = [manifest_path, result_path, evaluation_path, cond_lost_path]
        verify_input_fingerprints(data)
        output_root = Path(args.out).resolve()
        atomic_json(final_marker, dict(run_hash=run_hash, status="complete",
                    artifacts={str(p.resolve().relative_to(output_root)): sha_file(p) for p in artifacts}))
        primary = staar_result["u17"]["primary_bbj_designated"]
        log(f"STAAR flat ratio={primary['ratio']:.6f} z={primary['z']:+.4f}")
        log("RESULT " + json.dumps(dict(version="v9", arm="staar",
                                         evaluation_file=str(evaluation_path))))
        print("L1_DONE", flush=True)
        lock.close()
        return

    if args.kappa_eb_guide:
        trainer = Trainer(args, data, None, residuals, run_hash, truth=truth)
        frozen = trainer.freeze_head(inner_tr)
        guide = kappa_eb_guide(trainer, inner_tr, frozen, list(args.kappa_grid))
        guide_path = run_dir / "kappa_eb_guide.json"
        atomic_json(guide_path, guide)
        artifacts = [manifest_path, guide_path, cond_lost_path]
        verify_input_fingerprints(data)
        output_root = Path(args.out).resolve()
        atomic_json(final_marker, dict(run_hash=run_hash, status="complete",
                    artifacts={str(p.resolve().relative_to(output_root)): sha_file(p) for p in artifacts}))
        log(f"KAPPA_EB_GUIDE khat={guide['kappa_hat']:.6f} curvature={guide['curvature_at_hat']:.6g} fits=0")
        print("L1_DONE kappa_eb_guide", flush=True)
        lock.close()
        return

    design = Design(args, data, inner_tr)
    trainer = Trainer(args, data, design, residuals, run_hash, truth=truth)
    inner_frozen = trainer.freeze_head(inner_tr)
    trainer.gi, trainer.frozen = inner_tr, inner_frozen
    gradcheck = trainer.gradcheck(inner_tr, inner_frozen) if args.smoke else None
    if args.gradcheck_only:
        atomic_json(Path(args.out) / (tag + ".gradcheck.done"), gradcheck)
        print("L1_DONE gradcheck", flush=True)
        lock.close()
        return

    kappa_selection = None
    if args.v9_enabled:
        kappa_eval = Evaluation(trainer, inner_va, inner_frozen, args.inner_null_perm,
                                args.inner_seed, known, truth=truth)
        kappa_candidates, kappa_null, kappa_lam = {}, {}, {}
        for kappa in args.kappa_grid:
            trainer.kappa = float(kappa)
            trainer.gi, trainer.frozen = inner_tr, inner_frozen
            _, gradient = trainer.data_loss(design.initial)
            strong = design.lambda_strong(gradient)
            require(np.isfinite(strong), "nonfinite kappa lambda strong")
            if strong <= ZERO_SD:
                a_k, opt = design.initial.copy(), None
            else:
                a_k, opt = trainer.fit(strong, inner_tr, inner_frozen,
                                       args.inner_maxiter, tag + f"_kappa{kappa:g}_strong")
            f_k = design.f_of(a_k)
            phi_k = np.exp(np.clip(f_k, -CL, CL)).astype(np.float32)
            candidate, kappa_null[float(kappa)] = kappa_eval.evaluate(phi_k)
            candidate.update(kappa=float(kappa), lam_strong=float(strong), optimizer=opt,
                             rms_f=float(np.sqrt(np.square(f_k[design.train_pairs]).mean())))
            kappa_candidates[float(kappa)] = candidate
            kappa_lam[float(kappa)] = float(strong)
            metric = candidate["u17"]["primary_bbj_designated"]
            log(f"kappa={kappa:g} lambda_strong={strong:.6g} U17 ratio={metric['ratio']:.6f} z={metric['z']:+.4f}")
        chosen_kappa, kappa_selection = select_kappa_one_se(kappa_candidates, kappa_null)
        for key, row in kappa_selection["candidates"].items():
            kappa = float(key)
            metric = kappa_candidates[kappa]["u17"]["primary_bbj_designated"]
            row.update(lam_strong=kappa_lam[kappa], ratio=metric["ratio"],
                       null_sd=metric["null_sd"])
        kappa_selection.update(null_plan=kappa_eval.plan,
                               objective="v8 phenotype prior-only" if args.arm_prior_only else "BBJ continuous target",
                               train_genes=int(len(inner_tr)), heldout_genes=int(len(inner_va)))
        trainer.kappa = chosen_kappa
        lam_strong = kappa_lam[chosen_kappa]
        kappa_selection_path = run_dir / "kappa_selection.json"
        atomic_json(kappa_selection_path, kappa_selection)
        del kappa_eval
    else:
        chosen_kappa = 0.0
        _, gradient = trainer.data_loss(design.initial)
        lam_strong = design.lambda_strong(gradient)
        require(np.isfinite(lam_strong), "nonfinite lambda strong")
        kappa_selection_path = None

    evaluation = Evaluation(trainer, inner_va, inner_frozen, args.inner_null_perm,
                            args.inner_seed, known, truth=truth)
    candidates, null_z = {}, {}
    candidates["inf"], null_z["inf"] = evaluation.evaluate()
    candidates["inf"].update(lam=None, divisor=None, rms_f=0.0, optimizer=None)
    if args.learning_curve:
        require(all(0 < f <= 1 for f in args.learning_curve), "learning-curve fractions in (0,1]")
        rng = np.random.default_rng(args.inner_seed + 7)
        perm = rng.permutation(len(inner_tr))
        lc = dict(fold=args.fold, arm=args.arm, fractions=list(args.learning_curve), lam_strong=float(lam_strong),
                  inner_tr_full=int(len(inner_tr)), inner_va=int(len(inner_va)), inner_maxiter=args.inner_maxiter,
                  fixed=dict(design="basis/knots/standardization from full inner_tr", head="pi/tau2 frozen on full inner_tr",
                             null="same inner_seed permutations as the full run -> z comparable to selection.json"),
                  subsample="nested: first round(frac*n) of one fixed permutation", rows=[])
        for frac in args.learning_curve:
            k = int(round(frac * len(inner_tr))); require(k >= 50, "subsample too small")
            gi = np.sort(np.asarray(inner_tr)[perm[:k]])
            for divisor in list(args.lam_divisors):
                lam = lam_strong / divisor; key = "div" + f"{divisor:g}"
                a_k, opt = trainer.fit(lam, gi, inner_frozen, args.inner_maxiter, tag + f"_lc{frac:g}_" + key)
                f = design.f_of(a_k); phi = np.exp(np.clip(f, -CL, CL)).astype(np.float32)
                cand, _ = evaluation.evaluate(phi)
                rs = cand["u17"]["primary_bbj_designated"] if args.v9_enabled else cand["ratio_stats"]["mean"]
                observed_ratio = rs["ratio"] if args.v9_enabled else rs["obs"]
                lc["rows"].append(dict(frac=frac, n_genes=int(k), key=key, lam=float(lam), obs=observed_ratio, null_mean=rs["null_mean"],
                                       null_sd=rs["null_sd"], z=rs["z"], rms_f=float(np.sqrt(np.square(f[design.train_pairs]).mean())),
                                       nit=(opt or {}).get("nit"), amax=float(np.abs(a_k).max()),
                                       per_trait_ratio=(rs["per_trait"] if args.v9_enabled else
                                                        {t: cand["ratio_stats"][t]["obs"] for t in cand["ratio_stats"] if t != "mean"})))
                log(f"LC frac={frac:g} genes={k} {key}: z={rs['z']:+.4f} obs={observed_ratio:.6f}")
                atomic_json(Path(args.out) / (tag + ".learning_curve.partial.json"), lc)
        atomic_json(Path(args.out) / (tag + ".learning_curve.json"), lc)
        print("L1_DONE learning_curve", flush=True)
        lock.close()
        return
    grid = list(args.lam_divisors) if args.lam is None else []
    if lam_strong <= ZERO_SD and args.lam is None:
        grid = []
    if args.lam is not None and math.isfinite(args.lam):
        grid = ["fixed"]
    pos = 0
    extension = False
    while pos < len(grid):
        divisor = grid[pos]
        lam = args.lam if divisor == "fixed" else lam_strong / divisor
        key = "fixed" if divisor == "fixed" else "div" + (f"{divisor:g}")
        a_k, opt = trainer.fit(lam, inner_tr, inner_frozen, args.inner_maxiter, tag + "_inner_" + key)
        f = design.f_of(a_k)
        phi = np.exp(np.clip(f, -CL, CL)).astype(np.float32)
        candidates[key], null_z[key] = evaluation.evaluate(phi)
        candidates[key].update(lam=lam, divisor=divisor, optimizer=opt,
                               rms_f=float(np.sqrt(np.square(f[design.train_pairs]).mean())),
                               amax=float(np.abs(a_k).max()))
        chosen_metric = (candidates[key]["u17"]["primary_bbj_designated"] if args.v9_enabled
                         else candidates[key]["ratio_stats"]["mean"])
        log(f"inner {key}: ratio z={chosen_metric['z']:.4f}")
        selector = select_one_se_u17 if args.v9_enabled else select_one_se
        _, interim_selection = selector(candidates, null_z)
        if divisor == 1000 and not args.custom_grid and interim_selection["best"] == key:
            grid.append(10000)
            extension = True
        atomic_json(Path(args.out) / (tag + ".selection.partial.json"),
                    dict(candidates=candidates, selection=interim_selection, null_plan=evaluation.plan))
        pos += 1
    selector = select_one_se_u17 if args.v9_enabled else select_one_se
    chosen_key, selection = selector(candidates, null_z)
    if args.lam is not None:
        chosen_key = "inf" if math.isinf(args.lam) else "fixed"
        selection.update(chosen=chosen_key, reason="explicit fixed-lambda sensitivity; nested selection bypassed",
                         flat_selected=chosen_key == "inf", fixed_lambda_override=True)
    chosen = candidates[chosen_key]["lam"]
    selection.update(edge_extended=extension,
                     edge_rule="extend once to /10000 if /1000 is the strict strongest-first maximum ratio z before 1-SE",
                     zero_gradient_flat_only=lam_strong <= ZERO_SD and args.lam is None,
                     null_plan=evaluation.plan, frozen_head=inner_frozen,
                     train_traits=residuals.names, holdout_traits_excluded=args.holdout_trait,
                     design_sha256=design.signature, inner_cluster=design.cluster,
                     design_fit="inner training annotations only")
    if args.v9_enabled:
        selection.update(kappa=chosen_kappa, kappa_selection=kappa_selection)
    atomic_json(Path(args.out) / (tag + ".selection.json"), dict(selection=selection, candidates=candidates))
    del evaluation, trainer, design
    if args.device == "cuda":
        torch.cuda.empty_cache()
    design = Design(args, data, train)
    trainer = Trainer(args, data, design, residuals, run_hash, truth=truth)
    trainer.kappa = chosen_kappa
    frozen = trainer.freeze_head(train)
    if chosen is None:
        a, optimizer = design.initial.copy(), None
    else:
        a, optimizer = trainer.fit(chosen, train, frozen, args.maxiter, tag + "_final")
    f = design.f_of(a)
    phi = np.exp(np.clip(f, -CL, CL)).astype(np.float32)
    require(chosen is not None or np.array_equal(phi, np.ones_like(phi)), "flat model must be exactly flat")
    if args.holdout_trait:
        del trainer
        all_residuals = Residuals(args, data, rows)
        trainer = Trainer(args, data, design, all_residuals, run_hash, truth=truth)
        trainer.kappa = chosen_kappa
        eval_frozen = trainer.freeze_head(train, existing=frozen)
    else:
        all_residuals = residuals
        eval_frozen = frozen
    evaluation = Evaluation(trainer, test, eval_frozen, args.null_perm, args.test_seed,
                            known, truth=truth)
    test_result, _ = evaluation.evaluate(None if chosen is None else phi)
    clamp = float((np.abs(f) >= CL).mean())
    test_pairs = np.flatnonzero(np.isin(data.gid, test))
    clamp_test = float((np.abs(f[test_pairs]) >= CL).mean())
    correlations = {}
    correlations_test = {}
    for name in DIAG:
        x = np.nan_to_num(data.col(name))
        correlations[name] = float(np.corrcoef(x, np.log(phi))[0, 1]) if x.std() > 0 and phi.std() > 0 else 0.0
        xt, ft = x[test_pairs], np.log(phi[test_pairs])
        correlations_test[name] = float(np.corrcoef(xt, ft)[0, 1]) if xt.std() > 0 and ft.std() > 0 else 0.0
    out = dict(version=version, code_sha256=code_sha, run_hash=run_hash,
               arm="prior_only" if args.arm_prior_only else args.arm,
               mode="smoke" if args.smoke else f"fold{args.fold}", chr=data.use_chr,
               pairs=len(data.keys), genes=data.ng, train_genes=len(train), inner_tr=len(inner_tr),
               inner_va=len(inner_va), test_genes=len(test), n_param=design.nparam, blocks=design.meta,
               pair_weight="1/n_genes", penalty=dict(form="lambda a'Pa; spline G+Omega/q, binary/indicator ridge, linear E[x^2], nn I",
               CL=CL, EPS_F=EPS_F, lam_strong=lam_strong, lam_chosen=chosen, inner=candidates),
               selection=selection, frozen_head=frozen, evaluation_frozen_head=eval_frozen,
               perm_r=args.perm_r, resid_kind=all_residuals.kind,
               traits=names, train_traits=residuals.names, holdout_traits=args.holdout_trait,
               residual_diagnostics=all_residuals.diagnostics(),
               residual_sha256=all_residuals.hashes,
               objective=("v8 phenotype loss with BBJ-mixed gene prior" if args.arm_prior_only else
                          "BBJ continuous target" if args.v9_enabled else args.obj), obj_cap=args.obj_cap,
               objective_definition=("sum_trait,gene pi_g * (abs(T_g)-target_g,trait)^2" if
                                     args.v9_enabled and not args.arm_prior_only else
                                     "mean across training genes and training traits; uncapped J always reported"),
               cap_definition="-min(log m,c); caps positive log-BF reward, not positive loss",
               nit=optimizer["nit"] if optimizer else 0, converged=optimizer["success"] if optimizer else True,
               loss=optimizer["fun"] if optimizer else None, optimizer=optimizer,
               phi_q=np.quantile(phi, [0, .01, .5, .99, 1]).tolist(), clamp_frac=clamp,
               clamp_frac_test=clamp_test, corr_logphi_diag=correlations,
               corr_logphi_diag_test=correlations_test, gradcheck=gradcheck,
               resource_guard=guard, sec=round(time.time() - T0),
               inferential_status="fixed-phi diagnostic; reused fold0 is exploratory; no full-procedure null claim")
    if args.v9_enabled:
        out.update(kappa=chosen_kappa, kappa_selection=kappa_selection,
                   bbj=truth.manifest(), cadd_join=data.cadd_join,
                   conditional=cond_diagnostics)
    out.update(test_result)
    out["verdicts"] = verdicts(args, test_result, clamp_test, correlations_test)
    if args.arm == "cluster":
        effects = np.asarray(design.cluster["Q"]) @ a
        out["cluster"] = dict(design.cluster, effects=effects.tolist())
        atomic_npy(Path(args.out) / (tag + ".cluster_effects.npy"), effects)
    base = Path(args.out) / tag
    atomic_npy(Path(str(base) + ".a.npy"), a)
    atomic_npy(Path(str(base) + ".phi.npy"), phi)
    export_basis(Path(str(base) + ".basis.pkl"), design, code_sha)
    atomic_json(Path(str(base) + ".coef_names.json"), dict(coef_names=design.names))
    atomic_json(Path(str(base) + ".json"), out)
    evaluation_path = None
    if args.v9_enabled:
        evaluation_path = run_dir / "evaluation.json"
        atomic_json(evaluation_path, dict(run_hash=run_hash, fold=args.fold, smoke=args.smoke,
                                          primary_decision_only=True, result=test_result,
                                          u17=test_result["u17"]))
    artifacts = [Path(str(base) + suffix) for suffix in
                 (".a.npy", ".phi.npy", ".basis.pkl", ".coef_names.json", ".json", ".selection.json", ".manifest.json")]
    if args.v9_enabled:
        artifacts.extend([evaluation_path, kappa_selection_path, cond_lost_path])
    if args.arm == "cluster":
        artifacts.append(Path(str(base) + ".cluster_effects.npy"))
    verify_input_fingerprints(data)
    output_root = Path(args.out).resolve()
    atomic_json(Path(str(base) + ".done"), dict(run_hash=run_hash, status="complete",
                artifacts={str(p.resolve().relative_to(output_root)): sha_file(p) for p in artifacts}))
    log("RESULT " + json.dumps(dict(version=version, result_file=str(base) + ".json",
                                   chosen=chosen_key, ratio_mean=out["ratio_mean"])))
    print("L1_DONE", flush=True)
    lock.close()

if __name__ == "__main__":
    try:
        main()
    except RuntimeError as error:
        raise SystemExit(str(error)) from None
    except (ValueError, KeyError) as error:
        raise SystemExit("GATE FAIL: input parse/schema error (" + type(error).__name__ + ")") from None
