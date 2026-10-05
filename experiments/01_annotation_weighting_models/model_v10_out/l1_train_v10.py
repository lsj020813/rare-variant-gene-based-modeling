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


import argparse
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys
import time

OWNER = Path(_config_path('${PROJECT_ROOT}'))
WORKSPACE = Path(__file__).resolve().parent
READ_ROOTS = (WORKSPACE, OWNER / 'work/ref', OWNER / 'work/ref15')
FOLDS = {0: {'7', '13', '16', '19', '21'},
         1: {'2', '8', '9', '10', '15'}, 2: {'1', '6', '14'},
         3: {'3', '12', '17', '22'}, 4: {'4', '5', '11', '18', '20'}}
SEALED = {'21', '22'}
S1 = ['cons', 'epi_active', 'epi_repr', 'epi_trans', 'tf', 'linsight', 'gpn_msa', 'dist_tss']
S2 = ['re2g_max', 'gh_elem_score', 'gh_link_score']
S3 = ['is_cage_prom', 'in_body', 'in_tss3kb', 'in_re2g', 't1_na', 'is_indel']
STAAR_ANNOTATIONS = ['cadd', 'cons', 'epi_active', 'epi_repr', 'epi_trans', 'tf', 'linsight', 'gpn_msa']
BANDS = ('ge5', '1-5', '0.1-1', 'lt0.1')
PRIMARY_BAND = '0.1-1'
C2_COLUMNS = ['maf_bbj', 'rsq_bbj', 'region_n', 'marginal_abs_z_rank']
EPS_F, NKNOT, NUM_EPS, ZERO_SD = 0.05, 6, 1e-8, 1e-12
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

def atomic_json(path, value):
    path = owned_path(path)
    tmp = owned_path(path.with_name(path.name + ".tmp"))
    with open(tmp, "w") as fh:
        json.dump(value, fh, ensure_ascii=False, indent=2, allow_nan=False)
        fh.write("\n")
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)

def atomic_tsv(path, fieldnames, rows):
    path = owned_path(path)
    tmp = owned_path(path.with_name(path.name + ".tmp"))
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
    tmp = owned_path(path.with_name(path.name + ".tmp"))
    with open(tmp, "wb") as fh:
        np.savez(fh, **arrays)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)

def atomic_npy(path, value):
    path = owned_path(path)
    tmp = owned_path(Path(str(path) + ".tmp"))
    with open(tmp, "wb") as fh:
        np.save(fh, value, allow_pickle=False)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)

def resource_guard(args):
    sys.dont_write_bytecode = True
    import fcntl
    import resource
    import subprocess
    require(1 <= args.threads <= 4, "v10 threads must be 1..4")
    require(0 < args.memory_gb <= 64, "v10 memory cap must be <=64 GiB")
    tmp = owned_path(args.tmpdir)
    tmp.mkdir(parents=True, exist_ok=True)
    require(tmp.is_dir(), "TMPDIR missing")
    allowed = sorted(os.sched_getaffinity(0))
    os.sched_setaffinity(0, allowed[:args.threads])
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
    concurrent_heavy = set()
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
                        if any(b.startswith(heavy) for b in basenames):
                            if getattr(args, "allow_concurrent_heavy", False):
                                concurrent_heavy.add(pid)
                            else:
                                require(False, "another owner heavy stage is active")
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
    if concurrent_heavy:
        print(f"[guard] E6B_PATCH: proceeding alongside {len(concurrent_heavy)} owner heavy-stage process(es) (--allow-concurrent-heavy)", flush=True)
    return lock, dict(threads=args.threads, memory_cap_bytes=cap, tmpdir=str(tmp), concurrent_owner_heavy_processes=len(concurrent_heavy),
                      available_ram_bytes=available, other_users=len(users),
                      guard="v10 stdlib usage guard; CPU affinity and <=8 GiB RLIMIT_AS; external stage process detection")

def owned_path(path):
    p = Path(path).resolve()
    require(within(p, WORKSPACE), "write outside v10 workspace")
    parent = p
    while not parent.exists():
        parent = parent.parent
    require(parent.stat().st_uid == os.getuid(), "write ancestor owner mismatch")
    return p

def within(path, root):
    try:
        Path(path).relative_to(root)
        return True
    except ValueError:
        return False

def input_path(path):
    p = Path(path).resolve()
    require(any(within(p, r.resolve()) for r in READ_ROOTS), "input outside allowed annotation reference roots")
    require(p.is_file(), "required annotation/PIP input is absent")
    return p

def load_libraries(threads):
    global np, minimize, null_space, SplineTransformer, expit, logit, norm, rankdata, spearmanr
    import numpy as np
    from scipy.optimize import minimize
    from scipy.linalg import null_space
    from scipy.special import expit, logit
    from scipy.stats import norm, rankdata, spearmanr
    from sklearn.preprocessing import SplineTransformer
    from threadpoolctl import threadpool_limits
    globals()['_thread_limiter'] = threadpool_limits(limits=threads)

def atomic_gzip_tsv(path, fieldnames, rows):
    path = owned_path(path)
    tmp = owned_path(path.with_name(path.name + '.tmp'))
    with gzip.open(tmp, 'wt', newline='') as fh:
        writer = csv.DictWriter(fh, delimiter='\t', fieldnames=fieldnames, lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    with open(tmp, 'rb') as fh:
        os.fsync(fh.fileno())
    os.replace(tmp, path)

def tsv_rows(path, required):
    path = input_path(path)
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        fields = reader.fieldnames or []
        require(len(fields) == len(set(fields)) and set(required) <= set(fields), 'TSV schema/duplicate headers')
        for row in reader:
            require(None not in row and all(x is not None for x in row.values()), 'TSV row width')
            if 'chr' in row:
                require(chromosome(row['chr']) not in SEALED, 'sealed chr21/22 in input; producer must exclude these rows')
            yield row

def chromosome(value):
    s = str(value).strip().lower()
    s = s[3:] if s.startswith('chr') else s
    require(re.fullmatch(r'[1-9][0-9]*', s) is not None and 1 <= int(s) <= 22, 'invalid autosomal chromosome')
    return str(int(s))

def integer(value, label):
    require(re.fullmatch(r'-?[0-9]+', str(value).strip()) is not None, label + ' must use integer text')
    return int(value)

def number(value, label, missing=False):
    if missing and str(value).strip().lower() in ('', 'na', 'nan', '.', 'null'):
        return float('nan')
    try:
        v = float(value)
    except (ValueError, TypeError):
        require(False, label + ' is not numeric')
    require(math.isfinite(v), label + ' must be finite')
    return v

def variant_key(value, ch, pos):
    parts = str(value).strip().split(':')
    require(len(parts) == 4, 'variant key must be chr:pos:ref:alt')
    c = chromosome(parts[0])
    p = integer(parts[1], 'key position')
    alleles = [s.upper() for s in parts[2:]]
    require(c == ch and p == pos and p > 0, 'variant key/chr/position disagree')
    require(all(re.fullmatch(r'[ACGTN]+', a) for a in alleles) and alleles[0] != alleles[1], 'variant allele schema')
    return ':'.join([c, str(p)] + alleles)

def read_traits(path):
    rows = []
    for row in tsv_rows(path, ['trait', 'trait_group', 'N']):
        name, group = row['trait'].strip(), row['trait_group'].strip()
        require(re.fullmatch(r'[A-Za-z0-9_]+', name) is not None, 'unsafe trait code')
        require(bool(group), 'empty trait group')
        if str(row['N']).strip().upper() in ('NA', ''):
            n = None
        else:
            n = integer(row['N'], 'trait N')
            require(n > 0, 'trait N must be positive')
        rows.append(dict(trait=name, trait_group=group, N=n))
    require(rows and len({r['trait'] for r in rows}) == len(rows), 'empty/duplicate traits')
    return rows

def band_of(maf):
    require(np.isfinite(maf).all() and ((maf >= 0) & (maf <= .5)).all(), 'MAF outside [0,0.5]')
    return np.select([maf >= .05, maf >= .01, maf >= .001], BANDS[:3], default=BANDS[3])

class Data:
    def __init__(self, args):
        self.args = args
        self.fingerprints = {}
        self.traits = read_traits(args.traits)
        self.fingerprint(args.traits)
        stats_path = input_path(args.annot_stats)
        self.fingerprint(stats_path)
        with open(stats_path) as fh:
            cols = json.load(fh).get('cols')
        require(isinstance(cols, list) and len(set(cols)) == len(cols), 'annotation stats cols schema')
        if len(cols) == 33 and 'cadd' not in cols:
            cols = cols + ['cadd']
        require(len(cols) == 34 and set(S1 + S2 + S3 + ['cadd']) <= set(cols), 'expected 34 annotation columns')
        self.cols = cols
        self.X, self.keys, self.chrs, self.pos38, self.in_band, self.maf, self.rsq = [], [], [], [], [], [], []
        seen, duplicates = set(), 0
        required = ['variant_hg19', 'chr', 'pos_hg19', 'pos_hg38', 'in_our_band', 'maf_bbj', 'rsq_bbj'] + cols
        for row in tsv_rows(args.bbj_annot, required):
            ch = chromosome(row['chr'])
            pos = integer(row['pos_hg19'], 'annotation position')
            key = variant_key(row['variant_hg19'], ch, pos)
            duplicates += key in seen
            seen.add(key)
            x = [number(row[c], 'annotation', missing=True) for c in cols]
            in_band = integer(row['in_our_band'], 'in_our_band')
            require(in_band in (0, 1), 'in_our_band must be binary')
            p38 = integer(row['pos_hg38'], 'hg38 position')
            require(p38 > 0, 'invalid hg38 position')
            self.keys.append(key)
            self.X.append(x)
            self.chrs.append(ch)
            self.pos38.append(p38)
            self.in_band.append(in_band)
            self.maf.append(number(row['maf_bbj'], 'BBJ MAF'))
            self.rsq.append(number(row['rsq_bbj'], 'BBJ rsq'))
        require(self.keys and duplicates == 0, 'annotation variant key duplicates/empty input')
        self.fingerprint(args.bbj_annot)
        self.X = np.asarray(self.X, dtype=np.float64)
        self.keys, self.chrs = np.asarray(self.keys), np.asarray(self.chrs)
        self.pos38, self.in_band = np.asarray(self.pos38), np.asarray(self.in_band, bool)
        self.maf, self.rsq = np.asarray(self.maf), np.asarray(self.rsq)
        self.variant_bands = band_of(self.maf)
        require(((self.rsq >= 0) & (self.rsq <= 1)).all(), 'rsq outside [0,1]')
        for c in S3 + (['is_typed'] if 'is_typed' in cols else []):
            require(np.isin(self.col(c), [0, 1]).all(), 'binary annotation contract')
        na = self.col('t1_na').astype(bool)
        for col in S1[:7] + ['cadd']:
            values = self.col(col)
            require((~np.isfinite(values[na]) | (values[na] == 0)).all(), 't1_na marked value must be missing or zero')
            self.X[na, self.cols.index(col)] = np.nan
        require(np.array_equal(~np.isfinite(self.col('cadd')), na), 'CADD missingness must follow t1_na (use NA, not unmarked raw zero)')
        lookup = {k: i for i, k in enumerate(self.keys)}
        vi, ti, regions, cs, labels, zproxy = [], [], [], [], [], []
        audit = []
        region_chr = {}
        for t, trait in enumerate(self.traits):
            name = trait['trait']
            pip_path = input_path(Path(args.bbj_pip) / (name + '.tsv.gz'))
            local_seen, composite_seen, misses, dup_v, dup_c, total = set(), set(), 0, 0, 0, 0
            marginal = None
            if args.bbj_marginal:
                marginal = self.read_marginal(Path(args.bbj_marginal) / (name + '.tsv.gz'))
            required = ['variant_hg19', 'chr', 'pos', 'allele1', 'allele2', 'pip', 'cs_id', 'region']
            for row in tsv_rows(pip_path, required):
                ch, pos = chromosome(row['chr']), integer(row['pos'], 'PIP position')
                key = variant_key(row['variant_hg19'], ch, pos)
                require({row['allele1'].upper(), row['allele2'].upper()} == set(key.split(':')[2:]), 'PIP alleles disagree with variant key')
                region, cid = row['region'].strip(), integer(row['cs_id'], 'CS id')
                require(bool(region) and cid >= -1, 'region/CS schema')
                pip = number(row['pip'], 'PIP')
                require(0 <= pip <= 1, 'PIP outside [0,1]')
                total += 1
                dup_v += key in local_seen
                composite = (region, cid, key)
                dup_c += composite in composite_seen
                local_seen.add(key)
                composite_seen.add(composite)
                rkey = (t, region)
                require(region_chr.setdefault(rkey, ch) == ch, 'region spans chromosomes')
                if marginal is None:
                    require({'beta_marginal', 'se'} <= set(row), 'C2 needs --bbj-marginal or beta_marginal,se in each PIP file; no substitute')
                    zz = marginal_z(row)
                else:
                    require((key, region) in marginal, 'PIP to marginal join miss')
                    zz = marginal[(key, region)]
                if key not in lookup:
                    misses += 1
                    continue
                vi.append(lookup[key]); ti.append(t); regions.append(region); cs.append(cid); labels.append(pip); zproxy.append(zz)
            require(total > 0 and dup_v == 0 and dup_c == 0, 'empty PIP or nonunique variant/trait key; resolve multi-CS membership upstream')
            require(misses == 0, 'PIP to annotation join incomplete; no label imputation or silent intersection')
            audit.append(dict(trait=name, rows=total, duplicate_variant_trait=dup_v,
                              duplicate_variant_region_cs_trait=dup_c, unmatched_annotation=misses,
                              join_key='canonical hg19 chr:integer-pos:uppercase-ref:alt; unique per trait',
                              matched=total - misses))
            self.fingerprint(pip_path)
        self.vi, self.ti = np.asarray(vi, np.int32), np.asarray(ti, np.int32)
        self.regions, self.cs = np.asarray(regions), np.asarray(cs, np.int32)
        self.y = np.asarray(labels, np.float64)
        self.chr = self.chrs[self.vi]
        self.bands = self.variant_bands[self.vi]
        self.groups = np.asarray([self.traits[i]['trait_group'] for i in self.ti])
        self.region_blocks = make_blocks([(int(t), r) for t, r in zip(self.ti, self.regions)])
        self.cs_blocks = make_blocks([(int(t), r, int(c)) for t, r, c in zip(self.ti, self.regions, self.cs)], self.cs != -1)
        region_n, rank = np.empty(len(self.y)), np.empty(len(self.y))
        zproxy = np.asarray(zproxy)
        for idx in self.region_blocks:
            region_n[idx] = len(idx)
            rank[idx] = rankdata(zproxy[idx], method='average') / len(idx)
        self.C2 = np.column_stack([self.maf[self.vi], self.rsq[self.vi], region_n, rank])
        self.audit = dict(annotation_rows=len(self.keys), duplicate_annotation_keys=duplicates,
                          annotation_unlabelled=int(len(self.keys) - len(set(self.vi))),
                          labelled_pairs=len(self.y), joins=audit, sealed_chromosomes='rejected at ingest',
                          region_count_definition='all labelled variants in each trait/region before train-band filters',
                          c2_proxy='★ within-(trait,region) ECDF rank of abs(beta_marginal/se); proxy, not measured LD')
        for path in (args.oracle_spec, args.apply_annot):
            if path:
                self.fingerprint(path)

    def col(self, name):
        return self.X[:, self.cols.index(name)].copy()

    def fingerprint(self, path):
        p = input_path(path)
        st = p.stat()
        self.fingerprints[str(p)] = dict(size=st.st_size, mtime_ns=st.st_mtime_ns, sha256=sha_file(p))

    def read_marginal(self, path):
        result, duplicates = {}, 0
        for row in tsv_rows(path, ['variant_hg19', 'chr', 'pos', 'region', 'beta_marginal', 'se']):
            ch, pos = chromosome(row['chr']), integer(row['pos'], 'marginal position')
            key = (variant_key(row['variant_hg19'], ch, pos), row['region'].strip())
            duplicates += key in result
            result[key] = marginal_z(row)
        require(result and duplicates == 0, 'empty/nonunique marginal variant/region key')
        self.fingerprint(path)
        return result

def marginal_z(row):
    beta, se = number(row['beta_marginal'], 'marginal beta'), number(row['se'], 'marginal SE')
    require(se > 0, 'marginal SE must be positive')
    z = abs(beta / se)
    require(math.isfinite(z), 'nonfinite marginal z')
    return z

def make_blocks(keys, mask=None):
    groups = {}
    for i, key in enumerate(keys):
        if mask is None or mask[i]:
            groups.setdefault(key, []).append(i)
    return [np.asarray(v, dtype=np.int32) for v in groups.values()]

def weights(data, indices, y=None, base=None):
    indices = np.asarray(indices, dtype=int)
    y = data.y if y is None else y
    if base is None:
        base = weights_base(data, indices)
    base = np.asarray(base, dtype=float)
    high = y[indices] >= .1
    final = base.copy()
    masses = []
    for c in (False, True):
        mask = high == c
        mass = float(base[mask].sum())
        masses.append(mass)
        if mass:
            final[mask] /= mass
    return final, dict(base_mass=float(base.sum()), class_base_mass=masses,
                      class_final_mass=[float(final[~high].sum()), float(final[high].sum())],
                      scope='pooled across traits within this partition', proportional_constant=1)

def weights_base(data, indices):
    indices = np.asarray(indices, dtype=int)
    base = np.ones(len(indices), dtype=float)
    blocks = make_blocks([(int(data.ti[i]), data.regions[i], int(data.cs[i])) for i in indices], data.cs[indices] != -1)
    for block in blocks:
        base[block] = 1 / len(block)
    return base

def split_data(args, data):
    outer_chr = FOLDS[args.fold] - SEALED
    valid_chr = set().union(*(FOLDS[(args.fold + j) % 5] for j in range(1, args.inner_va_folds + 1))) - SEALED
    held = np.isin(data.groups, args.holdout_trait_group)
    require(set(args.holdout_trait_group) <= set(data.groups), 'unknown held-out trait group')
    outer = np.isin(data.chr, list(outer_chr))
    valid = np.isin(data.chr, list(valid_chr))
    eligible = np.ones(len(data.y), bool) if args.train_bands == 'all' else data.bands == args.train_bands
    train = np.flatnonzero(~held & ~outer & eligible)
    inner_tr = np.flatnonzero(~held & ~outer & ~valid & eligible)
    inner_va = np.flatnonzero(~held & ~outer & valid)
    test = np.flatnonzero(~held & outer)
    group_test = np.flatnonzero(held & ~outer)
    joint_test = np.flatnonzero(held & outer)
    require(all(len(x) for x in (train, inner_tr, inner_va, test)), 'empty chromosome/trait split')
    require(np.any(data.bands[inner_va] == PRIMARY_BAND), 'inner validation has no primary frequency band')
    train_t = set(data.ti[inner_tr])
    require(train_t == set(data.ti[train]) == set(data.ti[inner_va]) == set(data.ti[test]), 'seen traits must occur in every chromosome fitting/evaluation partition')
    if args.holdout_trait_group:
        require(len(group_test) and len(joint_test), 'empty trait-group or joint holdout')
    result = dict(train=train, inner_train=inner_tr, inner_valid=inner_va,
                  chromosome=test, trait_group=group_test, joint=joint_test)
    require(not set(inner_tr) & set(inner_va) and not set(train) & set(test), 'split overlap')
    return result

def phi_feature_set(cols, feature_mode, control):
    if control or feature_mode == 'all34':
        return list(cols)
    if feature_mode == 'v9-eighteen':
        return S1 + ['cadd'] + S2 + S3
    if feature_mode == 'v8-seventeen':
        return S1 + S2 + S3
    require(isinstance(feature_mode, str) and feature_mode.startswith('list:'), 'unknown --phi-columns mode')
    chosen = [c for c in feature_mode[5:].split(',') if c]
    require(chosen and len(set(chosen)) == len(chosen) and set(chosen) <= set(cols), 'explicit phi column list must be unique and within the 34 annotation columns')
    return chosen

class Design:
    def __init__(self, args, X, cols, train_rows, feature_mode='all34', control=False):
        self.args = args
        self.cols = list(cols)
        self.train_pairs = np.unique(train_rows)
        require(len(self.train_pairs) > 0, 'empty design fitting partition')
        self.is_train = np.zeros(len(X), bool)
        self.is_train[self.train_pairs] = True
        selected = phi_feature_set(self.cols, feature_mode, control)
        self.selected = list(selected)
        self.control = control
        self.basis, self.meta, self.blocks = {}, {}, []
        self.rules = {}
        binary_names = [] if control else [c for c in selected if c in S3 + ['is_typed']]
        binary = {c: X[:, self.cols.index(c)] for c in binary_names}
        for name in selected:
            x = np.asarray(X[:, self.cols.index(name)], dtype=float).copy()
            if name in binary_names:
                require(np.isin(x, [0, 1]).all(), 'invalid binary design column')
                self.binary(name, x)
                self.rules[name] = dict(source=name, transform='identity', indicator=False)
                continue
            mask = np.isfinite(x)
            transform = 'log1p' if name in ('dist_tss', 'gh_link_score') else 'identity'
            if transform == 'log1p':
                require((x[mask] >= 0).all(), 'negative log1p annotation')
                x[mask] = np.log1p(x[mask])
            self.continuous(name, x, mask)
            self.rules[name] = dict(source=name, transform=transform, indicator=False)
            need_indicator = name in S2 or (not control and name not in S1 + ['cadd'] and not mask.all())
            if need_indicator:
                dup = [c for c, b in binary.items() if np.array_equal(b[self.train_pairs].astype(bool), mask[self.train_pairs])]
                if dup:
                    self.meta[name]['indicator'] = 'duplicate of ' + dup[0]
                else:
                    self.binary(name + '_A', mask.astype(float), 'ind')
                    self.rules[name + '_A'] = dict(source=name, transform='identity', indicator=True)
                    self.meta[name]['indicator'] = name + '_A'
        self.B = np.hstack([b[2] for b in self.blocks])
        self.names = [f'{b[0]}:{j}' for b in self.blocks for j in range(b[2].shape[1])]
        self.order = [b[0] for b in self.blocks]
        self.P = np.zeros((len(self.names), len(self.names)))
        offset = 0
        for _, _, b, p in self.blocks:
            self.P[offset:offset + b.shape[1], offset:offset + b.shape[1]] = p
            offset += b.shape[1]
        self.blocks = None
        self.M = self.moment(self.B, self.train_pairs)
        self.nparam = self.B.shape[1]
        require(np.isfinite(self.B).all() and np.linalg.eigvalsh(self.P).min() > 0, 'invalid design/penalty')
        self.signature = digest(self.export())

    def moment(self, matrix, rows):
        result = np.zeros((matrix.shape[1], matrix.shape[1]))
        for lo in range(0, len(rows), self.args.design_chunk):
            block = matrix[rows[lo:lo + self.args.design_chunk]].astype(np.float64)
            result += block.T @ block
        return result / len(rows)

    def continuous(self, name, x, mask):
        tr = mask & self.is_train
        require(tr.any(), 'annotation column has no observed fitting rows')
        observed = x[tr]
        mu, sd = float(observed.mean()), float(observed.std()) + NUM_EPS
        xs = np.zeros(len(x))
        xs[mask] = (x[mask] - mu) / sd
        self.basis[name] = dict(mu=mu, sd=sd, kind='lin')
        meta = dict(n_obs_train=int(tr.sum()))
        self.meta[name] = meta
        knots = np.quantile(xs[tr], np.linspace(0, 1, NKNOT))
        unique = []
        for value in knots:
            if not unique or value - unique[-1] > 1e-6 * (knots[-1] - knots[0] + ZERO_SD):
                unique.append(value)
        if len(unique) < 3:
            size = float(np.square(xs[self.train_pairs]).mean())
            penalty = size if size > ZERO_SD else 1.0
            self.blocks.append((name, 'lin', xs[:, None].astype(np.float32), np.array([[penalty]])))
            meta.update(fallback='linear', zero_column_penalty=bool(size <= ZERO_SD))
            return
        knots = np.array(unique)
        spl = SplineTransformer(degree=3, knots=knots[:, None], include_bias=True,
                                extrapolation='constant').fit(xs[tr, None])
        raw = spl.transform(xs[:, None])
        mean_tr = raw[tr].mean(axis=0)
        Q = null_space(np.ones((1, raw.shape[1])))
        block = ((raw - mean_tr) @ Q) * mask[:, None]
        gram = block[tr].T @ block[tr] / tr.sum()
        require(np.linalg.matrix_rank(gram) == len(gram), 'v8 spline is unidentified on training values for column ' + name)
        bs = spl.bsplines_[0]
        grid = np.linspace(bs.t[3], bs.t[-4], 2001)
        deriv = bs.derivative(2)(grid)
        omega = Q.T @ (deriv.T @ deriv * (grid[1] - grid[0])) @ Q
        ev = np.linalg.eigvalsh(omega)
        rank = int((ev > NUM_EPS * ev.max()).sum())
        require(rank > 0, 'spline curvature rank zero')
        q = float(np.trace(np.linalg.solve(gram, omega)) / rank)
        require(np.isfinite(q) and q > 0, 'invalid curvature scaling')
        penalty = gram + omega / q
        require(np.linalg.eigvalsh(penalty).min() > 0, 'spline penalty not positive definite for column ' + name)
        self.blocks.append((name, 'spl', block.astype(np.float32), penalty))
        self.basis[name].update(kind='spl', knots=knots.tolist(), mean_tr=mean_tr.tolist(), Q=Q.tolist())
        meta.update(rank_omega=rank, q=q, n_basis=raw.shape[1])

    def binary(self, name, values, kind='bin'):
        p = float(values[self.train_pairs].mean())
        self.blocks.append((name, kind, (values - p)[:, None].astype(np.float32), np.array([[1.0]])))
        self.basis[name] = dict(kind=kind, p=p)

    def lambda_strong(self, gradient):
        pg = np.linalg.solve(self.P, gradient)
        return float(np.sqrt(max(0., pg @ self.M @ pg)) / (2 * EPS_F))

    def export(self):
        return dict(version='v10', cols=self.cols, selected_columns=self.selected,
                    inactive_columns=[c for c in self.cols if c not in self.selected],
                    basis=self.basis, meta=self.meta, order=self.order, rules=self.rules,
                    coef_names=self.names, NKNOT=NKNOT, degree=3, extrapolation='constant',
                    penalty='v8 G + Omega/q; Omega integrates squared second derivative',
                    missing='continuous centered basis zero; t1_na and applicability indicators as exported',
                    control=self.control, n_train_rows=len(self.train_pairs))

def transform_design(exported, X, cols):
    result = []
    for name in exported['order']:
        info, rule = exported['basis'][name], exported['rules'][name]
        x = np.asarray(X[:, cols.index(rule['source'])], dtype=float).copy()
        mask = np.isfinite(x)
        if rule['indicator']:
            x = mask.astype(float)
        if info['kind'] in ('bin', 'ind'):
            require(np.isfinite(x).all() and np.isin(x, [0, 1]).all(), 'invalid exported binary input')
            result.append((x - info['p'])[:, None])
            continue
        if rule['transform'] == 'log1p':
            require((x[mask] >= 0).all(), 'negative exported log1p input')
            x[mask] = np.log1p(x[mask])
        xs = np.zeros(len(x))
        xs[mask] = (x[mask] - info['mu']) / info['sd']
        if info['kind'] == 'lin':
            block = xs[:, None]
        else:
            knots = np.asarray(info['knots'])
            spl = SplineTransformer(degree=exported['degree'], knots=knots[:, None],
                                    include_bias=True, extrapolation=exported['extrapolation']).fit(knots[:, None])
            block = ((spl.transform(xs[:, None]) - np.asarray(info['mean_tr'])) @ np.asarray(info['Q'])) * mask[:, None]
        result.append(block)
    return np.hstack(result).astype(np.float32)

def _trim_memory():
    import gc, ctypes
    gc.collect()
    try:
        ctypes.CDLL('libc.so.6').malloc_trim(0)
    except Exception:
        pass

class Trainer:
    def __init__(self, args, data, design, train, control=False, y=None):
        self.args, self.data, self.design = args, data, design
        self.train, self.control = np.asarray(train), control
        self.y = data.y if y is None else y
        self.trait_ids = sorted(set(data.ti[self.train]))
        self.trait_offset = {t: i for i, t in enumerate(self.trait_ids)}
        self.group_names = sorted(set(data.groups[self.train])) if args.group_phi and not control else ['shared']
        self.q = design.nparam
        self.ncoef = self.q * len(self.group_names)
        self.P = np.kron(np.eye(len(self.group_names)), design.P)
        self.w, self.weight_meta = weights(data, self.train, self.y)
        if self.group_names == ['shared']:
            self.A = self.design.B[self.train if self.control else self.data.vi[self.train]]
        else:
            self.A = self.basis_for(self.train).astype(np.float32)
        self._step = max(1, args.design_chunk * 16)
        self.tidx = np.searchsorted(np.asarray(self.trait_ids), data.ti[self.train])
        require(np.array_equal(np.asarray(self.trait_ids)[self.tidx], data.ti[self.train]), 'trait index mapping')
        self.initial = np.zeros(self.ncoef + len(self.trait_ids))
        for t, j in self.trait_offset.items():
            mask = data.ti[self.train] == t
            mean = np.average(self.y[self.train][mask], weights=self.w[mask])
            self.initial[self.ncoef + j] = logit(np.clip(mean, NUM_EPS, 1 - NUM_EPS))

    def _Av(self, coef):
        out = np.empty(self.A.shape[0])
        for lo in range(0, self.A.shape[0], self._step):
            out[lo:lo + self._step] = self.A[lo:lo + self._step].astype(np.float64) @ coef
        return out

    def _Atv(self, v):
        out = np.zeros(self.A.shape[1])
        for lo in range(0, self.A.shape[0], self._step):
            out += self.A[lo:lo + self._step].astype(np.float64).T @ v[lo:lo + self._step]
        return out

    def basis_for(self, indices):
        base = self.design.B[indices if self.control else self.data.vi[indices]]
        if self.group_names == ['shared']:
            return base.astype(float)
        groups = self.data.groups[indices].copy()
        unknown = ~np.isin(groups, self.group_names)
        if unknown.any():
            require(self.args.unseen_phi_group in self.group_names, 'group-specific phi cannot predict unseen group without explicit --unseen-phi-group')
            groups = groups.astype(object)
            groups[unknown] = self.args.unseen_phi_group
        return np.hstack([base * (groups == g)[:, None] for g in self.group_names]).astype(float)

    def matrix(self, indices):
        basis = self.basis_for(indices)
        intercepts = np.column_stack([self.data.ti[indices] == t for t in self.trait_ids]).astype(float)
        return np.column_stack([basis, intercepts])

    def objective(self, a, lam):
        eta = self._Av(a[:self.ncoef]) + a[self.ncoef + self.tidx]
        loss = np.sum(self.w * (np.logaddexp(0., eta) - self.y[self.train] * eta))
        resid = self.w * (expit(eta) - self.y[self.train])
        grad = np.concatenate([self._Atv(resid), np.bincount(self.tidx, weights=resid, minlength=len(self.trait_ids))])
        pa = self.P @ a[:self.ncoef]
        loss += lam * (a[:self.ncoef] @ pa)
        grad[:self.ncoef] += 2 * lam * pa
        require(np.isfinite(loss) and np.isfinite(grad).all(), 'nonfinite logistic objective')
        return float(loss), grad

    def lambda_strong(self):
        _, g = self.objective(self.initial, 0.)
        if self.group_names == ['shared']:
            return self.design.lambda_strong(g[:self.ncoef])
        pg = np.linalg.solve(self.P, g[:self.ncoef])
        rms = np.sqrt(np.mean(np.square(self._Av(pg))))
        return float(rms / (2 * EPS_F))

    def fit(self, lam, maxiter):
        if lam is None or math.isinf(lam):
            return self.initial.copy(), dict(success=True, nit=0, intercept_only=True)
        p0 = expit(self._Av(self.initial[:self.ncoef]) + self.initial[self.ncoef + self.tidx])
        wpp = self.w * p0 * (1 - p0)
        q, n_t = self.ncoef, len(self.trait_ids)
        hessian = np.zeros((q + n_t, q + n_t))
        step = self._step
        for lo in range(0, len(wpp), step):
            Ac = self.A[lo:lo + step].astype(np.float64); WA = wpp[lo:lo + step, None] * Ac
            hessian[:q, :q] += Ac.T @ WA
            onehot = (self.tidx[lo:lo + step, None] == np.arange(n_t)[None, :]).astype(float)
            hessian[q:, :q] += onehot.T @ WA
        hessian[:q, q:] = hessian[q:, :q].T
        hessian[q:, q:] = np.diag(np.bincount(self.tidx, weights=wpp, minlength=n_t))
        hessian[:self.ncoef, :self.ncoef] += 2 * lam * self.P
        try:
            L = np.linalg.cholesky(hessian)
            transform = np.linalg.solve(L.T, np.eye(len(self.initial)))
            preconditioner = 'initial exact logistic Hessian Cholesky'
        except np.linalg.LinAlgError:
            transform = np.diag(1 / np.sqrt(np.maximum(np.diag(hessian), ZERO_SD)))
            preconditioner = 'initial Hessian diagonal; inherited numerical zero guard'
        def objective_scaled(z):
            loss, gradient = self.objective(self.initial + transform @ z, lam)
            return loss, transform.T @ gradient
        result = minimize(objective_scaled, np.zeros_like(self.initial), jac=True,
                          method='L-BFGS-B', options=dict(maxiter=maxiter, maxfun=3 * maxiter,
                          ftol=1e-9, gtol=1e-10, maxcor=20))
        require(result.status in (0, 1) and np.isfinite(result.x).all(), 'logistic optimizer failure')
        if self.args.require_convergence:
            require(result.success, 'optimizer did not converge in declared budget')
        fitted = self.initial + transform @ result.x
        require(np.isfinite(fitted).all(), 'nonfinite unscaled parameters')
        return fitted, dict(success=bool(result.success), nit=int(result.nit),
                            nfev=int(result.nfev), status=int(result.status), objective=float(result.fun),
                            preconditioner=preconditioner, objective_unchanged=True)

    def predict(self, a, indices):
        eta = self.basis_for(indices) @ a[:self.ncoef]
        unknown = np.ones(len(indices), bool)
        for t, j in self.trait_offset.items():
            mask = self.data.ti[indices] == t
            eta[mask] += a[self.ncoef + j]
            unknown[mask] = False
        if unknown.any():
            eta[unknown] += self.unseen_intercept(a)
        return open_probability(expit(eta))

    def unseen_intercept(self, a):
        if self.args.unseen_intercept is not None:
            return self.args.unseen_intercept
        require(self.args.unseen_intercept_rule == 'mean-trained', '★ unseen-trait intercept not sealed')
        return float(np.mean(a[self.ncoef:]))

def open_probability(p):
    return np.clip(p, np.nextafter(0., 1.), np.nextafter(1., 0.))

def ce_values(p, y):
    p = open_probability(np.asarray(p))
    return -(y * np.log(p) + (1 - y) * np.log1p(-p))

def staar_pair_weights(data, name, reference=None):
    values = data.col(name)
    require(reference is not None and len(reference), '★ C1 ECDF reference not supplied')
    reference = np.unique(reference)
    ref = values[reference]
    observed_ref = np.sort(ref[np.isfinite(ref)])
    require(len(observed_ref), 'C1 annotation has no observed reference values')
    p = np.full(len(values), .5)
    mask = np.isfinite(values)
    left = np.searchsorted(observed_ref, values[mask], side='left')
    right = np.searchsorted(observed_ref, values[mask], side='right')
    p[mask] = np.where(right > left, (left + right + 1) / 2, right) / len(reference)
    return p, dict(M=len(reference), n_observed=len(observed_ref),
                   n_missing=int((~np.isfinite(ref)).sum()), tie_method='average', missing_weight=.5,
                   denominator='all unique reference variants', fitted_labels=0)

def cauchy_combine(probabilities):
    p = np.asarray(open_probability(probabilities), dtype=np.longdouble)
    cauchy = np.mean(np.tan((np.longdouble(.5) - p) * np.longdouble(np.pi)), axis=0)
    result = np.longdouble(.5) - np.arctan(cauchy) / np.longdouble(np.pi)
    return open_probability(np.asarray(result, dtype=float))

def acat_equivalent_z(scores):
    p = 2 * norm.sf(np.abs(scores))
    result = norm.isf(np.maximum(cauchy_combine(p) / 2, np.nextafter(0., 1.)))
    require(np.isfinite(result).all(), 'nonfinite ACAT equivalent z')
    return result

def c1_predictions(args, data, train, indices):
    require(args.c1_ecdf_reference in ('training', 'evaluation'), '★ C1 ECDF reference not sealed')
    reference = np.unique(data.vi[train if args.c1_ecdf_reference == 'training' else indices])
    columns = data.cols if args.c1_columns == 'all34' else STAAR_ANNOTATIONS
    require(args.c1_columns in ('all34', 'v9-eight'), '★ C1 annotation set not sealed')
    predictions, meta = {}, {}
    for col in columns:
        pred, meta[col] = staar_pair_weights(data, col, reference)
        predictions['C1:' + col] = open_probability(pred[data.vi[indices]])
    predictions['C1'] = cauchy_combine(np.stack(list(predictions.values())))
    return predictions, dict(reference=args.c1_ecdf_reference, annotations=meta,
                            combination='equal-weight Cauchy score of ECDF predictions; STAAR-O reference only')

def c1_calibrate(args, data, train, y):
    c1, _ = c1_predictions(args, data, train, train)
    x = logit(np.clip(c1['C1'], NUM_EPS, 1 - NUM_EPS))
    w, _ = weights(data, train, y)
    yt = y[train]
    def obj(theta):
        eta = theta[0] + theta[1] * x
        loss = np.sum(w * (np.logaddexp(0., eta) - yt * eta))
        r = w * (expit(eta) - yt)
        return float(loss), np.array([r.sum(), (r * x).sum()])
    start = np.array([logit(np.clip(np.average(yt, weights=w), NUM_EPS, 1 - NUM_EPS)), 0.])
    res = minimize(obj, start, jac=True, method='L-BFGS-B', options=dict(maxiter=500, ftol=1e-12, gtol=1e-10))
    require(np.isfinite(res.x).all(), 'C1 calibration failed')
    return dict(intercept=float(res.x[0]), slope=float(res.x[1]), success=bool(res.success), nit=int(res.nit), n_train=int(len(train)),
                fitted_on='outer training rows (same rows as final C2/phi fits)', model='expit(a + b*logit(C1_cauchy))')

def apply_c1_calibration(cal, c1_scores):
    return open_probability(expit(cal['intercept'] + cal['slope'] * logit(np.clip(c1_scores, NUM_EPS, 1 - NUM_EPS))))

class PermutationPlan:
    def __init__(self, data, indices, count, seed):
        self.indices = np.asarray(indices)
        self.count, self.seed = count, seed
        self.blocks = make_blocks([(int(data.ti[i]), data.regions[i]) for i in indices])
        self.signature = digest(dict(seed=seed, count=count, indices_sha=hashlib.sha256(self.indices.tobytes()).hexdigest(),
                                     blocks=[len(b) for b in self.blocks], mode='within_trait_region'))

    def labels(self, y):
        rng = np.random.default_rng(self.seed)
        original = y[self.indices]
        for _ in range(self.count):
            result = original.copy()
            for block in self.blocks:
                result[block] = original[rng.permutation(block)]
            yield result

    def export(self):
        return dict(count=self.count, seed=self.seed, sha256=self.signature,
                    grouping='trait,region', same_for_all_arms=True,
                    recompute_class_weights=True, fixed_predictions=True,
                    interpretation='conditional held-out label null; no claim of refitted pipeline type-I calibration')

def null_summary(observed, null, direction=1):
    arr = np.asarray(null, float)
    good = arr[np.isfinite(arr)]
    if observed is None or not np.isfinite(observed) or len(good) < 2:
        return dict(observed=observed, null_mean=None, null_sd=None, z=None, p_greater=None,
                    n_valid=len(good), status='undefined_or_insufficient_null')
    mean, sd = float(good.mean()), float(good.std(ddof=1))
    return dict(observed=float(observed), null_mean=mean, null_sd=sd,
                z=float(direction * (observed - mean) / sd) if sd > ZERO_SD else None,
                p_greater=float((1 + np.sum(direction * good >= direction * observed)) / (len(good) + 1)),
                n_valid=len(good), status='ok' if sd > ZERO_SD else 'degenerate_null',
                null_mean_mcse=sd / math.sqrt(len(good)),
                null_values=[float(v) if np.isfinite(v) else None for v in arr])

def cs_spearman(predictions, labels, blocks):
    scores = [[] for _ in predictions]
    valid_cs = 0
    for b in blocks:
        if len(b) < 3:
            continue
        valid_cs += 1
        ry = rankdata(labels[b], method='average')
        ry -= ry.mean()
        sy = np.linalg.norm(ry)
        if sy <= ZERO_SD:
            continue
        for j, p in enumerate(predictions):
            rp = rankdata(p[b], method='average')
            rp -= rp.mean()
            sp = np.linalg.norm(rp)
            if sp > ZERO_SD:
                scores[j].append(float((rp @ ry) / (sp * sy)))
    return np.array([np.mean(s) if s else np.nan for s in scores]), [len(s) for s in scores], valid_cs

def calibration(p, y, w):
    cuts = np.quantile(p, np.linspace(0, 1, 11))
    assignment = np.searchsorted(cuts[1:-1], p, side='right')
    result = []
    for j in range(10):
        mask = assignment == j
        result.append(dict(decile=j + 1, n=int(mask.sum()), lower=float(cuts[j]), upper=float(cuts[j + 1]),
                           mean_p=float(p[mask].mean()) if mask.any() else None,
                           mean_pip=float(y[mask].mean()) if mask.any() else None,
                           weighted_mean_pip=float(np.average(y[mask], weights=w[mask])) if mask.any() else None))
    return result

def evaluate(args, data, indices, predictions, seed, y=None, only_bands=None):
    y = data.y if y is None else y
    indices = np.asarray(indices)
    names = list(predictions)
    require('C2' in names and 'flat' in names, 'evaluation requires C2 and flat')
    pred = np.stack([predictions[n] for n in names])
    require(pred.shape[1] == len(indices) and np.isfinite(pred).all(), 'prediction alignment')
    plan = PermutationPlan(data, indices, args.null_perm, seed)
    def label_iter():
        yield y[indices]
        yield from plan.labels(y)
    reports, raw = {}, {}
    bands = only_bands if only_bands is not None else BANDS
    for band in bands:
        mask = data.bands[indices] == band
        idx = indices[mask]
        if not len(idx):
            reports[band] = dict(status='empty_band', n=0, primary=band == PRIMARY_BAND)
            continue
        ps = pred[:, mask]
        blocks = make_blocks([(int(data.ti[i]), data.regions[i], int(data.cs[i])) for i in idx], data.cs[idx] != -1)
        ces, ranks, counts = [], [], []
        base_w = weights_base(data, idx)
        work_y = y.copy()
        for yy in label_iter():
            work_y[indices] = yy
            w, wm = weights(data, idx, work_y, base_w)
            ces.append(np.sum(ce_values(ps, yy[mask]) * w, axis=1) / w.sum())
            ranks_i, count_i, total = cs_spearman(ps, yy[mask], blocks)
            ranks.append(ranks_i)
            counts.append(count_i)
        ces, ranks = np.asarray(ces), np.asarray(ranks)
        delta = ces[:, [names.index('C2')]] - ces
        models = {}
        observed_w, observed_wm = weights(data, idx, y, base_w)
        for j, name in enumerate(names):
            rho = float(ranks[0, j]) if np.isfinite(ranks[0, j]) else None
            models[name] = dict(weighted_ce=null_summary(float(ces[0, j]), ces[1:, j], direction=-1),
                                delta_ce_c2=null_summary(float(delta[0, j]), delta[1:, j]),
                                cs_spearman=null_summary(rho, ranks[1:, j]),
                                cs_defined=counts[0][j], cs_eligible=total,
                                calibration=calibration(ps[j], y[idx], observed_w))
        ladder = {}
        sequence = [n for n in ('flat', 'C2', 'C1', 'phi') if n in names]
        if 'C1_cal' in names and 'C1' in names and 'phi' in names:
            sequence = sequence + [('C1', 'C1_cal'), ('C1_cal', 'phi')]
        pairs = [p for p in zip(sequence, sequence[1:]) if isinstance(p[0], str) and isinstance(p[1], str)] + [p for p in sequence if isinstance(p, tuple)]
        for left, right in pairs:
            diff = ces[:, names.index(left)] - ces[:, names.index(right)]
            ladder[left + '->' + right] = null_summary(float(diff[0]), diff[1:])
        reports[band] = dict(status='ok', n=len(idx), primary=band == PRIMARY_BAND,
                            weighting=observed_wm, models=models, ladder=ladder)
        raw[band] = dict(names=names, ce=ces, delta=delta, rank=ranks)
    return dict(bands=reports, null_plan=plan.export(),
                metrics='weighted mean soft-label CE; positive DeltaCE=CE(C2)-CE(arm); mean within-CS Spearman',
                weight_partition='B2 recomputed in each reported frequency band; shuffled labels reclassify weights',
                cs_undefined='constant predictions/PIP or CS with <3 band members excluded and counted'), raw

def select_one_se(candidates, null_scores):
    keys = list(candidates)
    best = max(keys, key=lambda k: candidates[k]['score'])
    best_score = candidates[best]['score']
    eligible, bounds = [], {}
    for key in keys:
        difference = np.asarray(null_scores[best]) - np.asarray(null_scores[key])
        se = float(difference.std(ddof=1))
        gap = best_score - candidates[key]['score']
        accept = bool(gap <= se + ZERO_SD)
        if accept:
            eligible.append(key)
        bounds[key] = dict(gap_to_best=float(gap), paired_null_se=se, eligible=accept)
    require(eligible, '1-SE found no candidate')
    return eligible[0], dict(best=best, chosen=eligible[0], strongest_first=keys, boundaries=bounds,
                            primary_band=PRIMARY_BAND, primary='raw held-out DeltaCE',
                            se_formula='sd_b(DeltaCE_best_null-DeltaCE_candidate_null), ddof=1; no sqrt(B)',
                            interpretation='inherited v8 conditional-null tolerance heuristic, not classical CV standard error')

def make_design(args, data, train, control=False):
    return Design(args, data.C2 if control else data.X,
                  C2_COLUMNS if control else data.cols,
                  train if control else np.unique(data.vi[train]), args.phi_columns, control)

def fit_select(args, data, design, train, valid, control, reference, out, tag, y=None, fixed_strong=None):
    y = data.y if y is None else y
    trainer = Trainer(args, data, design, train, control, y)
    strong = trainer.lambda_strong() if fixed_strong is None else fixed_strong
    grid = [('inf', None)] + [('div' + f'{d:g}', strong / d) for d in args.lam_divisors]
    if args.lam is not None:
        grid = [('fixed', args.lam)]
    candidates, nulls, models = {}, {}, {}
    flat = np.full(len(valid), np.mean(y[train]))
    for key, lam in grid:
        a, opt = trainer.fit(lam, args.inner_maxiter)
        pred = trainer.predict(a, valid)
        pp = dict(flat=flat, C2=flat if control else reference, candidate=pred)
        report, raw = evaluate(args, data, valid, pp, args.inner_seed, y, only_bands=[PRIMARY_BAND])
        r = raw[PRIMARY_BAND]
        j = r['names'].index('candidate')
        score = r['delta'][:, j]
        candidates[key] = dict(score=float(score[0]), lam=lam, optimizer=opt,
                               primary_metrics=report['bands'][PRIMARY_BAND]['models']['candidate'])
        nulls[key], models[key] = score[1:], a
        chosen, selection = select_one_se(candidates, nulls)
        partial = dict(candidates=candidates, selection=selection, design_sha256=design.signature,
                       lambda_strong=strong, null_plan=report['null_plan'],
                       c2_selection_reference='flat' if control else 'selected inner C2')
        if (key == 'div1000' and selection['best'] == key and args.lam is None and
                not getattr(args, 'custom_grid', False) and args.lam_divisors == [10., 100., 1000.]):
            grid.append(('div10000', strong / 10000.))
        partial['edge_extension'] = any(k == 'div10000' for k, _ in grid)
        atomic_json(Path(out) / (tag + '.selection.partial.json'), partial)
        log(f'{tag} {key}: inner fit complete, optimizer_success={opt["success"]}')
    atomic_json(Path(out) / (tag + '.selection.json'), partial)
    return dict(trainer=trainer, a=models[chosen], lam=candidates[chosen]['lam'],
                key=chosen, selection=partial, strong=strong)

def export_model(args, data, trainer, a):
    offsets = {data.traits[t]['trait']: float(a[trainer.ncoef + j]) for t, j in trainer.trait_offset.items()}
    global_intercept = float(np.mean(list(offsets.values())))
    default_intercept = None
    if args.export_intercept is not None:
        default_intercept = args.export_intercept
    elif args.export_trait is not None:
        require(args.export_trait in offsets, 'export trait was not trained')
        default_intercept = offsets[args.export_trait]
    elif args.export_intercept_rule == 'mean-trained':
        default_intercept = global_intercept
    return dict(version='v10', design=trainer.design.export(),
                coefficients={g: a[j * trainer.q:(j + 1) * trainer.q].tolist() for j, g in enumerate(trainer.group_names)},
                global_intercept=global_intercept, trait_intercepts=offsets,
                trait_offsets={k: v - global_intercept for k, v in offsets.items()},
                trait_groups={r['trait']: r['trait_group'] for r in data.traits if r['trait'] in offsets},
                default_intercept=default_intercept, default_phi_group=args.export_phi_group,
                output='sigmoid(additive annotation score + selected intercept), no exp/clamp burden transform',
                default_intercept_rule=args.export_intercept_rule, export_trait=args.export_trait,
                missing_convention=trainer.design.export()['missing'], code_sha256=sha_file(__file__),
                interpretation='class-balanced PIP regression score; sampling weights can alter marginal calibration')

class PhiPredictor:
    def __init__(self, model):
        if isinstance(model, (str, Path)):
            with open(input_path(model)) as fh:
                model = json.load(fh)
        require(model.get('version') == 'v10', 'export version')
        self.model = model

    def predict(self, X, cols, trait=None, intercept=None, phi_group=None):
        m = self.model
        X = np.asarray(X, dtype=float).copy()
        if 't1_na' in cols:
            na = X[:, cols.index('t1_na')] == 1
            for col in S1[:7] + ['cadd']:
                if col in cols:
                    values = X[:, cols.index(col)]
                    require((~np.isfinite(values[na]) | (values[na] == 0)).all(), 'application t1_na value must be missing or zero')
                    X[na, cols.index(col)] = np.nan
        groups = list(m['coefficients'])
        if groups == ['shared']:
            group = 'shared'
        else:
            group = phi_group or m['trait_groups'].get(trait) or m['default_phi_group']
            require(group in groups, '★ exported group-specific phi needs an explicit trained group')
        if trait is not None:
            require(trait in m['trait_intercepts'], 'unknown exported trait; supply an explicit intercept without trait')
            b = m['trait_intercepts'][trait]
        else:
            b = m['default_intercept'] if intercept is None else intercept
        require(b is not None and math.isfinite(b), '★ export/application intercept not sealed')
        basis = transform_design(m['design'], X, cols)
        return open_probability(expit(basis @ np.asarray(m['coefficients'][group]) + b))

    def apply_phi(self, annot_row, trait=None, intercept=None, phi_group=None):
        cols = self.model['design']['cols']
        x = np.array([[number(annot_row[c], 'application annotation', missing=True) for c in cols]])
        return float(self.predict(x, cols, trait, intercept, phi_group)[0])

def apply_phi(annot_row, model, trait=None, intercept=None, phi_group=None):
    return PhiPredictor(model).apply_phi(annot_row, trait, intercept, phi_group)

def export_band_scores(args, data, exported, out):
    require(args.apply_annot, '★ --apply-annot with validated variant_hg38 and band is required for band score export')
    predictor = PhiPredictor(exported)
    cols = data.cols
    X, keys, bands = [], [], []
    seen = set()
    for row in tsv_rows(args.apply_annot, ['variant_hg38', 'chr', 'pos_hg38', 'band'] + cols):
        ch, pos = chromosome(row['chr']), integer(row['pos_hg38'], 'application hg38 position')
        key = variant_key(row['variant_hg38'], ch, pos)
        require(key not in seen and row['band'] in BANDS, 'duplicate application key or invalid band')
        seen.add(key)
        keys.append(key); bands.append(row['band'])
        X.append([number(row[c], 'application annotation', missing=True) for c in cols])
    require(X, 'empty application annotation table')
    X, bands = np.asarray(X), np.asarray(bands)
    scores = predictor.predict(X, cols)
    paths = {}
    for band in BANDS:
        selected = np.flatnonzero(bands == band)
        if not len(selected):
            continue
        path = Path(out) / ('phi_scores_' + band + '.tsv.gz')
        atomic_gzip_tsv(path, ['variant_hg38', 'p'],
                        (dict(variant_hg38=keys[i], p=float(scores[i])) for i in selected))
        paths[band] = dict(path=str(path), n=len(selected), sha256=sha_file(path))
    return paths

def verify_input_fingerprints(data):
    for path, previous in data.fingerprints.items():
        p = input_path(path)
        st = p.stat()
        require(st.st_size == previous['size'] and st.st_mtime_ns == previous['mtime_ns'] and
                sha_file(p) == previous['sha256'], 'input changed during run')

def run_pipeline(args, data, out, y=None, extras=True):
    out = owned_path(out)
    out.mkdir(parents=True, exist_ok=True)
    y = data.y if y is None else y
    splits = split_data(args, data)
    tr, va = splits['inner_train'], splits['inner_valid']
    inner_c2 = make_design(args, data, tr, True)
    selected_c2 = fit_select(args, data, inner_c2, tr, va, True, None, out, 'C2', y)
    c2_valid = selected_c2['trainer'].predict(selected_c2['a'], va)
    inner_phi = make_design(args, data, tr)
    selected_phi = fit_select(args, data, inner_phi, tr, va, False, c2_valid, out, 'phi', y)
    train = splits['train']
    final_c2 = Trainer(args, data, make_design(args, data, train, True), train, True, y)
    c2_a, c2_opt = final_c2.fit(selected_c2['lam'], args.maxiter)
    final_phi = Trainer(args, data, make_design(args, data, train), train, False, y)
    phi_a, phi_opt = final_phi.fit(selected_phi['lam'], args.maxiter)
    final_phi.comparator = (final_c2, c2_a)
    exported = export_model(args, data, final_phi, phi_a)
    flat = float(np.mean(y[train]))
    c1_cal = c1_calibrate(args, data, train, y)
    reports = {}
    for axis in ('chromosome', 'trait_group', 'joint'):
        idx = splits[axis]
        if not len(idx):
            reports[axis] = dict(status='not_requested')
            continue
        predictions, c1_meta = c1_predictions(args, data, train, idx)
        predictions = dict(flat=np.full(len(idx), flat), C2=final_c2.predict(c2_a, idx),
                           **predictions, C1_cal=apply_c1_calibration(c1_cal, predictions['C1']), phi=final_phi.predict(phi_a, idx))
        reports[axis], _ = evaluate(args, data, idx, predictions, args.test_seed, y)
        reports[axis]['C1_specification'] = c1_meta
        reports[axis]['C1_calibration'] = c1_cal
    result = dict(version='v10', run_kind='external_PIP_annotation_only', fold=args.fold, folds_used={str(k): sorted(v) for k, v in FOLDS.items()},
                  train_bands=args.train_bands, holdout_trait_group=args.holdout_trait_group,
                  group_phi=args.group_phi, primary_band=PRIMARY_BAND,
                  split_counts={k: len(v) for k, v in splits.items()},
                  split_hashes={k: hashlib.sha256(v.tobytes()).hexdigest() for k, v in splits.items()},
                  feature_contract=dict(phi_columns=args.phi_columns, selected=final_phi.design.selected,
                                        inactive=final_phi.design.export()['inactive_columns']),
                  selection=dict(C2=selected_c2['selection'], phi=selected_phi['selection']),
                  final_fit=dict(C2=c2_opt, phi=phi_opt, C2_lambda=selected_c2['lam'], phi_lambda=selected_phi['lam']),
                  evaluation=reports, flat_training_mean_pip=flat, join_audit=data.audit,
                  input_fingerprints=data.fingerprints, code_sha256=sha_file(__file__),
                  annotation_only=True, target='continuous external PIP; no missing-region zero labels')
    if extras and args.learning_curves:
        result['learning_curves'] = learning_curves(args, data, splits, inner_phi, inner_c2,
                                                   selected_phi['strong'], selected_c2['strong'], out, y)
    if extras and args.export_phi:
        filename = args.export_phi
        path = Path(filename) if Path(filename).is_absolute() else out / filename
        atomic_json(path, exported)
        result['export_phi'] = dict(path=str(path), sha256=sha_file(path))
        if args.apply_annot:
            result['band_scores'] = export_band_scores(args, data, exported, out)
        else:
            result['band_scores'] = dict(status='pending_apply_annotation_contract',
                                        required='hg38-oriented variant key and annotation-only table via --apply-annot')
    atomic_npz(out / 'fitted_parameters.npz', phi=phi_a, C2=c2_a)
    atomic_json(out / 'RESULT.partial.json', result)
    return result, exported, final_phi, phi_a

def learning_curves(args, data, splits, phi_design, c2_design, phi_strong, c2_strong, out, y):
    tr, va = splits['inner_train'], splits['inner_valid']
    require(args.train_bands == 'all', 'two-axis learning curves require --train-bands all')
    rng = np.random.default_rng(args.inner_seed + 7)
    blocks = make_blocks([(int(data.ti[i]), data.regions[i], int(data.cs[i])) for i in tr], data.cs[tr] != -1)
    strata = {}
    for b in blocks:
        indices = tr[b]
        signature = (int(data.ti[indices[0]]), tuple(int(np.sum(data.bands[indices] == band)) for band in BANDS))
        strata.setdefault(signature, []).append(indices)
    ordered = [[v[j] for j in rng.permutation(len(v))] for v in strata.values()]
    negative = tr[data.cs[tr] == -1]
    points = []
    for fraction in (.25, .5, .75, 1.):
        kept = [b for groups in ordered for b in groups[:int(round(fraction * len(groups)))]]
        indices = np.sort(np.concatenate([negative] + kept))
        points.append(('signal', f'{fraction:g}', indices, dict(fraction=fraction, selected_cs=len(kept), total_cs=len(blocks))))
    for bands in ((PRIMARY_BAND,), (PRIMARY_BAND, '1-5'), (PRIMARY_BAND, '1-5', 'ge5')):
        indices = tr[np.isin(data.bands[tr], bands)]
        points.append(('band', '+'.join(bands), indices, dict(bands=list(bands))))
    result = dict(fixed=dict(design='full inner training unique variants/pairs; knots/scalers unchanged',
                             validation='same inner-valid variants', permutations='same inner_seed and region blocks',
                             lambda_grid='full inner-train lambda_strong/divisors for both fitted arms'),
                  cs_sampling='nested prefixes within trait x exact CS band-count-vector strata; round as v8',
                  non_cs='cs_id=-1 rows remain fixed in signal-axis curves', rows=[])
    for i, (axis, point, indices, meta) in enumerate(points):
        point_dir = owned_path(Path(out) / ('curve_' + str(i)))
        point_dir.mkdir(parents=True, exist_ok=True)
        require(len(indices) and set(data.ti[indices]) == set(data.ti[tr]), 'learning-curve point lacks training traits')
        fit_c2 = fit_select(args, data, c2_design, indices, va, True, None, point_dir, 'C2', y, c2_strong)
        c2_pred = fit_c2['trainer'].predict(fit_c2['a'], va)
        fit_phi = fit_select(args, data, phi_design, indices, va, False, c2_pred, point_dir, 'phi', y, phi_strong)
        c1, _ = c1_predictions(args, data, tr, va)
        pred = dict(flat=np.full(len(va), np.mean(y[indices])), C2=c2_pred, **c1,
                    phi=fit_phi['trainer'].predict(fit_phi['a'], va))
        metrics, _ = evaluate(args, data, va, pred, args.inner_seed, y, only_bands=[PRIMARY_BAND])
        row = dict(axis=axis, point=point, n_train=len(indices), **meta,
                   trait_band_counts={data.traits[t]['trait']: {b: int(np.sum((data.ti[indices] == t) & (data.bands[indices] == b))) for b in BANDS} for t in sorted(set(data.ti[tr]))},
                   subset_sha256=hashlib.sha256(indices.tobytes()).hexdigest(),
                   phi_design_sha256=phi_design.signature, C2_design_sha256=c2_design.signature,
                   evaluation=metrics)
        result['rows'].append(row)
        atomic_json(Path(out) / 'learning_curve.partial.json', result)
        atomic_json(point_dir / 'POINT_DONE', dict(row_sha256=digest(row)))
        log(f'learning curve {i + 1}/{len(points)} complete')
    atomic_json(Path(out) / 'learning_curve.json', result)
    return result

def oracle_signal(data, spec):
    if spec.get('phi_model'):
        require('intercepts' in spec, 'oracle trait intercepts missing')
        with open(input_path(spec['phi_model'])) as fh:
            m = json.load(fh)
        require(m.get('version') == 'v10' and list(m['coefficients']) == ['shared'], 'phi_model oracle needs a shared-phi v10 export')
        basis = transform_design(m['design'], data.X, data.cols)
        signal = basis @ np.asarray(m['coefficients']['shared'])
        require(set(spec['intercepts']) >= {r['trait'] for r in data.traits}, 'oracle trait intercepts missing')
        return signal
    require(set(spec) >= {'terms', 'intercepts'} and spec['terms'], 'oracle specification schema')
    signal = np.zeros(len(data.keys))
    for term in spec['terms']:
        require(set(term) >= {'column', 'coefficient', 'center', 'scale', 'power'}, '★ all oracle term values must be explicit')
        require(term['column'] in data.cols and term['scale'] > 0 and term['power'] in (1, 2, 3), 'oracle term domain')
        x = data.col(term['column'])
        x = np.where(np.isfinite(x), (x - term['center']) / term['scale'], 0.)
        signal += term['coefficient'] * x ** term['power']
    require(set(spec['intercepts']) >= {r['trait'] for r in data.traits}, 'oracle trait intercepts missing')
    return signal

def generate_pip(data, spec, delta, rng=None):
    law = spec.get('label_law')
    require(law in ('deterministic_sigmoid', 'beta_mean_sigmoid'), '★ unsupported/unsealed oracle PIP law')
    signal = oracle_signal(data, spec)
    intercept = np.asarray([spec['intercepts'][data.traits[t]['trait']] for t in data.ti])
    mean = open_probability(expit(delta * signal[data.vi] + intercept))
    if law == 'beta_mean_sigmoid':
        require(rng is not None and spec.get('concentration', 0) > 0, '★ beta PIP law needs explicit concentration and RNG')
        strength = spec['concentration']
        mean = rng.beta(mean * strength, (1 - mean) * strength)
    return open_probability(mean), signal

def classify(args, result, simulation=None):
    required = dict(equivalence_margin=args.equivalence_margin,
                    power_target=args.power_target, minimum_relevant_delta=args.minimum_relevant_delta)
    missing = [k for k, v in required.items() if v is None]
    if args.alpha is None:
        return dict(category='판정 불능', reason='significance level not sealed', unsealed=['alpha'])
    main = result['evaluation']['chromosome']['bands'][PRIMARY_BAND]
    if main['status'] != 'ok':
        return dict(category='판정 불능', reason='primary-band holdout is absent')
    ladder = main['ladder']
    c1_phi = ladder['C1->phi']
    phi_c2 = main['models']['phi']['delta_ce_c2']
    significant = lambda s: s['status'] == 'ok' and s['observed'] > 0 and s['p_greater'] <= args.alpha
    if significant(c1_phi):
        return dict(category='성공', reason='phi exceeds C1 on primary-band CE at the sealed alpha')
    if missing or not simulation or not simulation.get('power_adequate', False):
        return dict(category='판정 불능', reason='unsealed criteria or D3 has not established adequate power',
                    unsealed=missing, power_evidence=bool(simulation and simulation.get('power_adequate')))
    close = abs(c1_phi['observed']) <= args.equivalence_margin
    flat_ce = main['models']['flat']['weighted_ce']['observed']
    all_flat = all(abs(m['weighted_ce']['observed'] - flat_ce) <= args.equivalence_margin
                   for name, m in main['models'].items() if name in ('C2', 'C1', 'phi'))
    if (close and significant(phi_c2)) or all_flat:
        return dict(category='의미 있는 실패', reason='sealed descriptive CE margin plus C2 improvement, or all arms within flat margin',
                    equivalence_note='margin is descriptive; no formal equivalence claim from randomization SD')
    return dict(category='판정 불능', reason='specified success/meaningful-failure conditions not met; no fourth category invented')

def binomial_interval(successes, n, alpha):
    from scipy.stats import beta
    require(n > 0 and 0 < alpha < 1, 'binomial interval settings')
    lo = 0. if successes == 0 else float(beta.ppf(alpha / 2, successes, n - successes + 1))
    hi = 1. if successes == n else float(beta.ppf(1 - alpha / 2, successes + 1, n - successes))
    return [lo, hi]

def simulate(args, data, out):
    require(args.oracle_spec and args.effect_grid is not None and args.simulation_reps is not None,
            '★ --oracle-spec, --effect-grid and --simulation-reps must be sealed for D3')
    require(args.alpha is not None and args.power_target is not None, '★ simulation alpha/power target missing')
    with open(input_path(args.oracle_spec)) as fh:
        spec = json.load(fh)
    result = dict(generator=spec, effect_grid=args.effect_grid, replicates=args.simulation_reps,
                  null=[], oracle=[], power=[], delta_star_candidate=None,
                  seal='Delta* is a candidate only; human sealing required',
                  structure='original annotation, trait/region and CS membership; no cohort inputs')
    all_indices = np.arange(len(data.y))
    for mode in ('null', 'oracle', 'power'):
        deltas = args.effect_grid if mode == 'power' else [None]
        for j, delta in enumerate(deltas):
            rows = []
            for r in range(args.simulation_reps):
                seeds = np.random.SeedSequence([args.test_seed, ('null', 'oracle', 'power').index(mode), j, r]).generate_state(3)
                if mode == 'null':
                    plan = PermutationPlan(data, all_indices, 1, int(seeds[0]))
                    yy = next(plan.labels(data.y))
                    truth_signal = None
                else:
                    yy, truth_signal = generate_pip(data, spec, 1. if mode == 'oracle' else delta,
                                                   np.random.default_rng(int(seeds[0])))
                subargs = argparse.Namespace(**vars(args))
                subargs.inner_seed, subargs.test_seed = int(seeds[1]), int(seeds[2])
                subout = Path(out) / 'simulation' / mode / str(j) / str(r)
                rr, model, trainer, coeff = run_pipeline(subargs, data, subout, yy, extras=False)
                metric = rr['evaluation']['chromosome']['bands'][PRIMARY_BAND]['ladder']['C1->phi']
                passed = metric['status'] == 'ok' and metric['observed'] > 0 and metric['p_greater'] <= args.alpha
                row = dict(replicate=r, delta=delta, metric=metric, significant=bool(passed))
                if truth_signal is not None:
                    idx = split_data(args, data)['chromosome']
                    predicted = trainer.basis_for(idx) @ coeff[:trainer.ncoef]
                    target = truth_signal[data.vi[idx]]
                    row['phi_recovery_correlation'] = float(np.corrcoef(predicted, target)[0, 1]) if np.std(predicted) > ZERO_SD else None
                    cp, _ = c1_predictions(args, data, split_data(args, data)['train'], idx)
                    comparator, comparator_a = trainer.comparator
                    mean_spec = dict(spec, label_law='deterministic_sigmoid')
                    oracle_mean, _ = generate_pip(data, mean_spec, 1. if mode == 'oracle' else delta)
                    oracle_pred = dict(flat=np.full(len(idx), np.mean(yy[split_data(args, data)['train']])),
                                       C2=comparator.predict(comparator_a, idx), **cp,
                                       oracle=oracle_mean[idx], pointwise_pip_ceiling=yy[idx])
                    ceiling, _ = evaluate(subargs, data, idx, oracle_pred, subargs.test_seed, yy, [PRIMARY_BAND])
                    row['oracle_ceiling'] = ceiling['bands'][PRIMARY_BAND]['models']['oracle']
                    row['pointwise_pip_ceiling'] = ceiling['bands'][PRIMARY_BAND]['models']['pointwise_pip_ceiling']
                rows.append(row)
                atomic_json(subout / 'RESULT.json', rr)
                atomic_json(subout / 'L1_DONE', dict(result_sha256=sha_file(subout / 'RESULT.json')))
                log(f'D3 {mode}, grid point {j + 1}/{len(deltas)}, replicate {r + 1}/{args.simulation_reps} complete')
                atomic_json(Path(out) / 'simulation.partial.json', dict(**result, current=dict(mode=mode, delta=delta, rows=rows)))
            successes = sum(row['significant'] for row in rows)
            interval = binomial_interval(successes, len(rows), args.alpha)
            entry = dict(delta=delta, rows=rows, rejection_rate=successes / len(rows), interval=interval)
            result[mode].append(entry)
    null_ci = result['null'][0]['interval']
    result['type1_check'] = dict(compatible_with_alpha=null_ci[0] <= args.alpha <= null_ci[1],
                                 interval=null_ci, alpha=args.alpha,
                                 conclusion='finite-replicate compatibility only; precision must be assessed by human')
    eligible = [p['delta'] for p in result['power'] if p['interval'][0] >= args.power_target]
    result['delta_star_candidate'] = min(eligible) if eligible else None
    result['power_adequate'] = bool(result['delta_star_candidate'] is not None and
                                   args.minimum_relevant_delta is not None and
                                   result['delta_star_candidate'] <= args.minimum_relevant_delta and
                                   result['type1_check']['compatible_with_alpha'])
    stochastic = spec['label_law'] == 'beta_mean_sigmoid'
    result['power_interpretation'] = ('independent synthetic PIP draws conditional on fixed annotations/region/CS structure'
                                     if stochastic else 'deterministic oracle: conditional probability across permutation seeds only')
    type1_adequate = args.type1_max is not None and null_ci[1] <= args.type1_max
    result['power_adequate'] = bool(result['power_adequate'] and stochastic and type1_adequate)
    result['type1_check']['sealed_max'] = args.type1_max
    result['type1_check']['adequate_precision'] = type1_adequate
    atomic_json(Path(out) / 'simulation.json', result)
    return result

def safety_args(parser):
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--allow-concurrent-heavy', action='store_true', help='★ E6B_PATCH: downgrade owner heavy-stage overlap from failure to a logged note')
    parser.add_argument('--memory-gb', type=float, default=8)
    parser.add_argument('--tmpdir', default=str(WORKSPACE / '.scratch'))
    parser.add_argument('--device', choices=['cpu'], default='cpu')

def parser():
    ap = argparse.ArgumentParser(description=__doc__)
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument('--fold', type=int, choices=range(5))
    mode.add_argument('--make-smoke-fixture', metavar='DIR')
    mode.add_argument('--smoke-suite', metavar='FIXTURE_DIR')
    mode.add_argument('--apply-phi', metavar='MODEL_JSON', help='guarded annotation-only application of an exported model')
    ap.add_argument('--out', default=str(WORKSPACE / 'v10_results'))
    ap.add_argument('--bbj-annot')
    ap.add_argument('--bbj-pip')
    ap.add_argument('--bbj-marginal', help='optional side directory; otherwise beta_marginal/se required in PIP files')
    ap.add_argument('--traits', help='TSV headers: trait,trait_group,N')
    ap.add_argument('--annot-stats', default=str(OWNER / 'work/ref15/annot/cache/fm_all.stats.json'))
    ap.add_argument('--phi-columns', help="★ 'all34' | 'v9-eighteen' (S1+cadd+S2+S3=18) | 'v8-seventeen' (S1+S2+S3, no cadd) | 'list:c1,c2,...' explicit; no default  [E6_PATCH 2026-09-08]")
    ap.add_argument('--group-phi', action='store_true', help='group-specific additive phi; default off')
    ap.add_argument('--holdout-trait-group', action='append', default=[])
    ap.add_argument('--unseen-intercept', type=float, help='★ explicit unseen-trait intercept')
    ap.add_argument('--unseen-intercept-rule', choices=['mean-trained'], help='★ alternative unseen-trait intercept rule')
    ap.add_argument('--unseen-phi-group', help='★ explicit trained group used to transfer group-specific phi')
    ap.add_argument('--train-bands', choices=['all', '1-5', '0.1-1'], default='all')
    ap.add_argument('--inner-va-folds', type=int, choices=[1, 2], default=2)
    ap.add_argument('--folds-override', help='E6G_PATCH: JSON {fold: [chromosomes]} replacing FOLDS (subset simulations); SEALED unchanged')
    ap.add_argument('--lam-divisors', type=float, nargs='+', default=[10., 100., 1000.])
    ap.add_argument('--lam', type=float, help='fixed-lambda sensitivity; inf means intercept-only')
    ap.add_argument('--null-perm', type=int, default=50)
    ap.add_argument('--inner-seed', type=int, default=20260908)
    ap.add_argument('--test-seed', type=int, default=123)
    ap.add_argument('--inner-maxiter', type=int, default=30)
    ap.add_argument('--maxiter', type=int, default=60)
    ap.add_argument('--design-chunk', type=int, default=16384)
    ap.add_argument('--require-convergence', action='store_true')
    ap.add_argument('--c1-columns', choices=['all34', 'v9-eight'], help='★ C1 annotation set; no default')
    ap.add_argument('--c1-ecdf-reference', choices=['training', 'evaluation'], help='★ ECDF reference population; no default')
    ap.add_argument('--learning-curves', action='store_true')
    ap.add_argument('--simulate', action='store_true')
    ap.add_argument('--oracle-spec', help='★ JSON: explicit terms/intercepts/PIP generative law')
    ap.add_argument('--effect-grid', type=float, nargs='+', help='★ delta grid; no default')
    ap.add_argument('--simulation-reps', type=int, help='★ independent simulation replicates; no default')
    ap.add_argument('--alpha', type=float, help='★ significance level; no default')
    ap.add_argument('--equivalence-margin', type=float, help='★ descriptive CE equivalence margin; no default')
    ap.add_argument('--power-target', type=float, help='★ minimum power; no default')
    ap.add_argument('--minimum-relevant-delta', type=float, help='★ effect at which power must be adequate')
    ap.add_argument('--type1-max', type=float, help='★ acceptable upper bound for null rejection rate')
    ap.add_argument('--export-phi', nargs='?', const='phi.json', help='export model JSON; optional relative filename')
    ap.add_argument('--export-trait', help='★ trained trait intercept to bind for transfer')
    ap.add_argument('--export-intercept', type=float, help='★ explicit transfer intercept')
    ap.add_argument('--export-intercept-rule', choices=['mean-trained'], help='★ transfer intercept rule')
    ap.add_argument('--export-phi-group', help='★ trained phi group for unlabelled transfer')
    ap.add_argument('--apply-annot', help='annotation-only TSV.gz: variant_hg38,chr,pos_hg38,band plus 34 annotations')
    safety_args(ap)
    return ap

def validate_args(args):
    require(args.null_perm >= 2, 'at least two permutations required for null SD')
    require(args.inner_maxiter > 0 and args.maxiter > 0 and args.design_chunk > 0, 'positive optimizer/chunk budgets')
    require(args.lam_divisors and all(math.isfinite(d) and d > 0 for d in args.lam_divisors) and
            args.lam_divisors == sorted(set(args.lam_divisors)), 'lambda divisors must be finite positive unique strongest-first')
    require(args.lam is None or args.lam >= 0, 'lambda must be nonnegative')
    for value in (args.alpha, args.power_target, args.type1_max):
        require(value is None or 0 < value < 1, 'probability decision setting outside (0,1)')
    for value in (args.equivalence_margin, args.minimum_relevant_delta):
        require(value is None or math.isfinite(value) and value >= 0, 'invalid nonnegative decision setting')
    require(args.effect_grid is None or args.effect_grid and all(math.isfinite(d) and d >= 0 for d in args.effect_grid) and
            args.effect_grid == sorted(set(args.effect_grid)), 'delta grid must be nonnegative unique ascending')
    require(args.simulation_reps is None or args.simulation_reps > 0, 'simulation replicate count')
    require(sum(v is not None for v in (args.export_trait, args.export_intercept, args.export_intercept_rule)) <= 1,
            'choose one export intercept convention')
    require(not (args.unseen_intercept is not None and args.unseen_intercept_rule is not None), 'choose one unseen intercept convention')
    for value in (args.export_intercept, args.unseen_intercept):
        require(value is None or math.isfinite(value), 'intercept must be finite')
    if args.apply_phi:
        require(args.apply_annot, '--apply-phi requires --apply-annot')
    if args.fold is not None:
        require(all([args.bbj_annot, args.bbj_pip, args.traits]), 'annotation/PIP/trait inputs required')
        require(all([args.phi_columns, args.c1_columns, args.c1_ecdf_reference]),
                '★ seal --phi-columns, --c1-columns, --c1-ecdf-reference before fitting')
        if args.holdout_trait_group:
            require(args.unseen_intercept is not None or args.unseen_intercept_rule is not None, '★ unseen-trait intercept convention required')
        if args.apply_annot:
            require(args.export_phi and any(v is not None for v in (args.export_trait, args.export_intercept, args.export_intercept_rule)),
                    '★ score export needs model export and bound transfer intercept')

def make_smoke_fixture(args):
    target = owned_path(args.make_smoke_fixture)
    require(not target.exists(), 'fixture target already exists')
    for sub in ('bbj_pip', 'bbj_marginal'):
        (target / sub).mkdir(parents=True, exist_ok=True)
    with open(input_path(args.annot_stats)) as fh:
        cols = json.load(fh)['cols']
    if 'cadd' not in cols:
        cols = cols + ['cadd']
    require(len(cols) == 34, 'fixture annotation schema')
    rng = np.random.default_rng(20260908)
    n = 20000
    X = rng.normal(size=(n, len(cols)))
    for name in S3 + ['is_typed']:
        X[:, cols.index(name)] = rng.integers(0, 2, n)
    for name in ('dist_tss', 'gh_link_score', 'n_genes', 'gh_n_genes'):
        X[:, cols.index(name)] = np.abs(X[:, cols.index(name)])
    missing = rng.random(n) < .05
    X[:, cols.index('t1_na')] = missing
    for name in S1[:7] + ['cadd']:
        X[missing, cols.index(name)] = np.nan
    for name in S2:
        X[rng.random(n) < .1, cols.index(name)] = np.nan
    chrs = (np.arange(n) // 1000 + 1).astype(str)
    pos = np.arange(n) + 100000
    keys = np.asarray([f'{c}:{p}:A:C' for c, p in zip(chrs, pos)])
    band_index = np.tile(np.repeat(np.arange(4), 250), 20)
    bounds = ((.05, .5), (.01, .05), (.001, .01), (0., .001))
    maf = np.empty(n)
    for b, (lower, upper) in enumerate(bounds):
        mask = band_index == b
        maf[mask] = rng.uniform(lower, upper, mask.sum())
    rsq = rng.uniform(.5, 1., n)
    traits = [dict(trait='SYN_A', trait_group='지질', N=10000),
              dict(trait='SYN_B', trait_group='당대사', N=10000),
              dict(trait='SYN_C', trait_group='혈압', N=10000)]
    spec = dict(label_law='beta_mean_sigmoid', concentration=1000.,
                terms=[dict(column='cons', coefficient=1.2, center=0., scale=1., power=1),
                       dict(column='cadd', coefficient=-1.0, center=0., scale=1., power=1),
                       dict(column='epi_active', coefficient=.5, center=0., scale=1., power=2)],
                intercepts=dict(SYN_A=-3., SYN_B=-3.4, SYN_C=-2.7))
    atomic_json(target / 'oracle_spec.json', spec)
    atomic_json(target / 'fm_all.stats.json', dict(cols=cols, synthetic=True))
    atomic_tsv(target / 'traits.tsv', ['trait', 'trait_group', 'N'], traits)
    def annotation_rows():
        for i in range(n):
            yield dict(variant_hg19=keys[i], chr=chrs[i], pos_hg19=int(pos[i]), pos_hg38=int(pos[i] + 1000000),
                       **{c: float(X[i, j]) if np.isfinite(X[i, j]) else 'NA' for j, c in enumerate(cols)},
                       in_our_band=int(.01 <= maf[i] < .05), maf_bbj=float(maf[i]), rsq_bbj=float(rsq[i]))
    atomic_gzip_tsv(target / 'bbj_annot.tsv.gz', ['variant_hg19', 'chr', 'pos_hg19', 'pos_hg38'] + cols + ['in_our_band', 'maf_bbj', 'rsq_bbj'], annotation_rows())
    class Synthetic:
        def col(self, name):
            return self.X[:, self.cols.index(name)]
    synthetic = Synthetic()
    synthetic.X, synthetic.cols, synthetic.keys, synthetic.traits = X, cols, keys, traits
    synthetic.vi, synthetic.ti = np.tile(np.arange(n), 3), np.repeat(np.arange(3), n)
    pip, truth = generate_pip(synthetic, spec, 1., rng)
    pip_fields = ['variant_hg19', 'chr', 'pos', 'allele1', 'allele2', 'pip', 'cs_id', 'region']
    marginal_fields = ['variant_hg19', 'chr', 'pos', 'region', 'beta_marginal', 'se']
    for t, trait in enumerate(traits):
        pip_rows, marginal_rows = [], []
        for i in range(n):
            local, r = i % 1000, i // 1000
            if local >= 900:
                continue
            cid = -1
            if t == r % 3:
                if local < 100:
                    cid = 0
                elif 500 <= local < 600:
                    cid = 1
            region = f'SYN_REGION_{r}'
            pip_rows.append(dict(variant_hg19=keys[i], chr=chrs[i], pos=int(pos[i]), allele1='A', allele2='C',
                                 pip=float(pip[t * n + i]), cs_id=cid, region=region))
            marginal_rows.append(dict(variant_hg19=keys[i], chr=chrs[i], pos=int(pos[i]), region=region,
                                      beta_marginal=float(rng.normal()), se=1.))
        atomic_gzip_tsv(target / 'bbj_pip' / (trait['trait'] + '.tsv.gz'), pip_fields, pip_rows)
        atomic_gzip_tsv(target / 'bbj_marginal' / (trait['trait'] + '.tsv.gz'), marginal_fields, marginal_rows)
    def application_rows():
        for i in range(n):
            if .01 <= maf[i] < .05:
                yield dict(variant_hg38=f'{chrs[i]}:{pos[i] + 1000000}:A:C', chr=chrs[i],
                           pos_hg38=int(pos[i] + 1000000), band='1-5',
                           **{c: float(X[i, j]) if np.isfinite(X[i, j]) else 'NA' for j, c in enumerate(cols)})
    atomic_gzip_tsv(target / 'apply_annot.tsv.gz', ['variant_hg38', 'chr', 'pos_hg38', 'band'] + cols, application_rows())
    atomic_npy(target / 'oracle_signal.npy', truth)
    fixture_files = sorted(p for p in target.rglob('*') if p.is_file())
    manifest = dict(synthetic=True, variants=n, traits=3, regions=20, cs=40,
                    files={str(p.relative_to(target)): sha_file(p) for p in fixture_files},
                    smoke_settings=dict(fold=0, phi_columns='all34', c1_columns='v9-eight', c1_ecdf_reference='training',
                                        holdout_trait_group=['혈압'], unseen_intercept_rule='mean-trained',
                                        export_intercept_rule='mean-trained', null_perm=50,
                                        inner_maxiter=120, maxiter=240, learning_curves=True,
                                        alpha=.05, equivalence_margin=.01, power_target=.8,
                                        minimum_relevant_delta=1., effect_grid=[0., 1.], simulation_reps=2),
                    smoke_gates=dict(oracle_correlation_min=.9, null_z_mean_tolerance=.5, null_replicates=20),
                    settings_status='★ synthetic engineering choices only; NOT real-data defaults/seals')
    atomic_json(target / 'FIXTURE.json', manifest)
    log('synthetic fixture created: 20000 variants, 3 traits, 20 regions, 40 CS')
    return manifest

def fixture_args(args, target):
    target = owned_path(target)
    with open(input_path(target / 'FIXTURE.json')) as fh:
        manifest = json.load(fh)
    require(manifest.get('synthetic') is True, 'smoke requires synthetic fixture manifest')
    for relative, expected in manifest['files'].items():
        path = input_path(target / relative)
        require(within(path, target) and sha_file(path) == expected, 'fixture checksum/provenance mismatch')
    cfg = argparse.Namespace(**vars(args))
    for name, value in manifest['smoke_settings'].items():
        setattr(cfg, name, value)
    cfg.bbj_annot, cfg.bbj_pip = str(target / 'bbj_annot.tsv.gz'), str(target / 'bbj_pip')
    cfg.bbj_marginal, cfg.traits = str(target / 'bbj_marginal'), str(target / 'traits.tsv')
    cfg.annot_stats = str(target / 'fm_all.stats.json')
    cfg.oracle_spec, cfg.apply_annot = str(target / 'oracle_spec.json'), str(target / 'apply_annot.tsv.gz')
    cfg.export_phi, cfg.smoke_suite = 'phi.json', None
    validate_args(cfg)
    return cfg, manifest

def smoke_suite(args, guard):
    from v10_smoke import run_smoke
    return run_smoke(sys.modules[__name__], args, guard)

def apply_exported(args, guard):
    from types import SimpleNamespace
    out = owned_path(args.out)
    require(not (out / 'L1_DONE').exists(), 'application output is already completed; use a new directory')
    paths = [input_path(args.apply_phi), input_path(args.apply_annot)]
    fingerprints = {str(p): sha_file(p) for p in paths}
    with open(paths[0]) as fh:
        model = json.load(fh)
    require(model.get('version') == 'v10', 'application model version')
    if args.export_intercept is not None:
        model['default_intercept'] = args.export_intercept
    elif args.export_trait is not None:
        require(args.export_trait in model['trait_intercepts'], 'application trait absent from export')
        model['default_intercept'] = model['trait_intercepts'][args.export_trait]
        model['default_phi_group'] = model['trait_groups'][args.export_trait]
    elif args.export_intercept_rule == 'mean-trained':
        model['default_intercept'] = model['global_intercept']
    if args.export_phi_group:
        model['default_phi_group'] = args.export_phi_group
    scores = export_band_scores(args, SimpleNamespace(cols=model['design']['cols']), model, out)
    require(all(sha_file(p) == fingerprints[str(p)] for p in paths), 'application inputs changed during run')
    result = dict(version='v10', kind='annotation_only_application', scores=scores,
                  guard=guard, input_sha256=fingerprints, code_sha256=sha_file(__file__))
    atomic_json(out / 'RESULT.json', result)
    atomic_json(out / 'L1_DONE', dict(result_sha256=sha_file(out / 'RESULT.json')))
    print('L1_DONE apply_phi', flush=True)

def main():
    args = parser().parse_args()
    args.custom_grid = '--lam-divisors' in sys.argv
    if getattr(args, 'folds_override', None):
        override = json.loads(args.folds_override)
        require(sorted(int(k) for k in override) == [0, 1, 2, 3, 4], 'folds-override must define folds 0..4')
        FOLDS.clear(); FOLDS.update({int(k): set(str(c) for c in v) for k, v in override.items()})
        print(f'[E6G_PATCH] FOLDS overridden: {sorted((k, sorted(v)) for k, v in FOLDS.items())}', flush=True)
    validate_args(args)
    lock, guard = resource_guard(args)
    try:
        load_libraries(args.threads)
        if args.apply_phi:
            apply_exported(args, guard)
            return
        if args.make_smoke_fixture:
            manifest = make_smoke_fixture(args)
            atomic_json(Path(args.make_smoke_fixture) / 'L1_DONE', dict(kind='fixture', manifest_sha256=digest(manifest), guard=guard))
            print('L1_DONE fixture', flush=True)
            return
        if args.smoke_suite:
            smoke_suite(args, guard)
            return
        out = owned_path(args.out)
        config = dict(arguments=vars(args), code_sha256=sha_file(__file__))
        done = out / 'L1_DONE'
        if done.exists():
            with open(done) as fh:
                saved = json.load(fh)
            require(saved.get('config_sha256') == digest(config), 'completed run configuration differs; use a new output directory')
            require(saved['result_sha256'] == sha_file(out / 'RESULT.json'), 'completed RESULT hash mismatch')
            data = Data(args); _trim_memory()
            with open(out / 'RESULT.json') as fh:
                old = json.load(fh)
            require(old['input_fingerprints'] == data.fingerprints, 'completed run input fingerprints changed')
            print('L1_DONE validated completed run', flush=True)
            return
        if (out / 'RUN.json').exists():
            with open(out / 'RUN.json') as fh:
                prior = json.load(fh)
            require(prior['arguments'] == config['arguments'] and prior['code_sha256'] == config['code_sha256'],
                    'incomplete output has a different configuration/code; use a new directory')
        atomic_json(out / 'RUN.json', dict(**config, guard=guard))
        data = Data(args); _trim_memory()
        result, _, _, _ = run_pipeline(args, data, out)
        simulation = simulate(args, data, out) if args.simulate else None
        if simulation:
            result['simulation'] = simulation
        result['verdict'] = classify(args, result, simulation)
        result['guard'] = guard
        verify_input_fingerprints(data)
        atomic_json(out / 'RESULT.json', result)
        atomic_json(done, dict(config_sha256=digest(config), result_sha256=sha_file(out / 'RESULT.json')))
        print('L1_DONE', flush=True)
    except Exception as error:
        message = str(error) if isinstance(error, RuntimeError) and str(error).startswith('GATE FAIL:') else type(error).__name__
        log(message)
        atomic_json(Path(args.out) / 'FAILED.json', dict(error=message, code_sha256=sha_file(__file__)))
        raise SystemExit(1) from None
    finally:
        lock.close()

if __name__ == '__main__':
    main()
