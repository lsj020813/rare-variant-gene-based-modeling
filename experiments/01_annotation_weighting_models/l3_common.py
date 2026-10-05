#!/usr/bin/env python3
import os as _cfg_os
import math as _cfg_math

def _config_number(name, cast=int, positive=False):
    value = cast(_cfg_os.environ[name])
    if not _cfg_math.isfinite(value) or (positive and value <= 0):
        raise ValueError(name + " must be finite" + (" and positive" if positive else ""))
    return value

import os as _os
N_SAMPLES = _config_number("N_SAMPLES", int, True)
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time

SEALED = frozenset()
S1 = ['cons', 'epi_active', 'epi_repr', 'epi_trans', 'tf', 'linsight', 'gpn_msa', 'dist_tss']
S2 = ['re2g_max', 'gh_elem_score', 'gh_link_score']
S3 = ['is_cage_prom', 'in_body', 'in_tss3kb', 'in_re2g', 't1_na', 'is_indel']
STAAR_ANNOTATIONS = ['cadd', 'cons', 'epi_active', 'epi_repr', 'epi_trans', 'tf', 'linsight', 'gpn_msa']
PHI_V9_EIGHTEEN = S1 + ['cadd'] + S2 + S3
NUM_EPS = 1e-8
CACHE_NCOL = 33
COHORT_N = N_SAMPLES
T0 = time.time()

np = sp = None

def require(condition, message):
    if not condition:
        raise RuntimeError('GATE FAIL: ' + message)

def log(message):
    print(f'[{time.time() - T0:8.1f}s] {message}', flush=True)

def load_libraries(threads):
    global np, sp
    import numpy
    import scipy.sparse
    from threadpoolctl import threadpool_limits
    np, sp = numpy, scipy.sparse
    globals()['_thread_limiter'] = threadpool_limits(limits=threads)

def sha_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()

def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, default=str).encode()).hexdigest()

def atomic_json(path, value):
    path = Path(path)
    tmp = path.with_name(path.name + '.tmp')
    with open(tmp, 'w') as fh:
        json.dump(value, fh, indent=1, sort_keys=True, default=str)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)

def atomic_tsv(path, header, rows):
    path = Path(path)
    tmp = path.with_name(path.name + '.tmp')
    with open(tmp, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(header)
        w.writerows(rows)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)

def norm_iid(value):
    s = str(value).strip()
    return s.split('_')[-1] if '_' in s else s

def import_v10(script):
    script = Path(script).resolve()
    require(script.is_file(), 'v10 trainer script not found')
    sys.path.insert(0, str(script.parent))
    mod = __import__(script.stem)
    require(hasattr(mod, 'PhiPredictor') and hasattr(mod, 'transform_design'),
            'v10 trainer lacks PhiPredictor/transform_design')
    mod.SEALED = set()
    return mod

class BandCache:

    def __init__(self, root, cadd_dir, cadd_source, chroms, t1_partial_policy):
        self.root = Path(root)
        self.fingerprints = {}
        path = self.root / 'annot/cache/fm_all.npz'
        self._fp(path)
        with np.load(path, allow_pickle=True) as z:
            X = z['X']
            self.cols = list(z['cols'].astype(str))
            keys = z['key37'].astype(str)
            genes = z['gene'].astype(str)
            chrs = z['chr'].astype(int).astype(str)
        require(X.shape[1] == CACHE_NCOL, 'band cache column count')
        require(len(set(self.cols)) == len(self.cols), 'duplicate cache columns')
        stats = self.root / 'annot/cache/fm_all.stats.json'
        self._fp(stats)
        with open(stats) as fh:
            require(list(json.load(fh)['cols']) == self.cols, 'cache/stats column disagreement')
        summaries = sorted((self.root / 'groupfiles_bwg').glob('chr*.summary.json'))
        require(bool(summaries), 'group summary files absent')
        expected = sum(int(json.load(open(p))['pairs']) for p in summaries)
        require(X.shape[0] == expected, 'cache rows vs groupfile pair total')
        self.n_pairs_all = int(X.shape[0])
        self.chroms_all = sorted(set(chrs), key=int)
        require(not (set(self.chroms_all) & set(SEALED)), 'sealed chromosome present (SEALED must be empty)')

        cadd = np.full(len(keys), np.nan, np.float32)
        self.cadd_chroms = []
        self.cadd_source = cadd_source
        if cadd_dir is not None:
            require(cadd_source in ('features', 'extract'), 'unknown --cadd-source')
            for ch in self.chroms_all:
                table = {}
                if cadd_source == 'features':
                    f = Path(cadd_dir) / f'chr{ch}.features.tsv.gz'
                    if not (f.is_file() and Path(cadd_dir, f'chr{ch}.done').is_file()):
                        continue
                    self._fp(f)
                    import gzip
                    with gzip.open(f, 'rt') as fh:
                        hdr = fh.readline().rstrip('\n').split('\t')
                        require(hdr[0] == 'key37' and 'cadd_phred' in hdr, 'features header')
                        j = hdr.index('cadd_phred')
                        for line in fh:
                            a = line.rstrip('\n').split('\t')
                            table[a[0]] = a[j]
                else:
                    f = Path(cadd_dir) / f'chr{ch}.cadd.tsv'
                    if not (f.is_file() and Path(str(f) + '.done').is_file()):
                        continue
                    self._fp(f)
                    with open(f) as fh:
                        for line in fh:
                            a = line.rstrip('\n').split('\t')
                            if len(a) >= 2 and a[0] != 'key37':
                                table[a[0]] = a[1]
                sel = np.flatnonzero(chrs == ch)
                vals = [table.get(k) for k in keys[sel]]
                cadd[sel] = [np.nan if v in (None, '', 'NA', 'nan', '.') else float(v) for v in vals]
                self.cadd_chroms.append(ch)
        self.cadd_available = sorted(self.cadd_chroms, key=int)

        X = np.hstack([X, cadd[:, None]])
        self.cols = self.cols + ['cadd']
        require(len(self.cols) == CACHE_NCOL + 1, 'design column count after cadd merge')

        na = X[:, self.cols.index('t1_na')] == 1
        partial = 0
        for cname in S1[:7] + ['cadd']:
            j = self.cols.index(cname)
            bad = na & np.isfinite(X[:, j]) & (X[:, j] != 0)
            partial += int(bad.sum())
            if t1_partial_policy == 'blank':
                X[bad, j] = np.nan
        require(t1_partial_policy == 'blank' or partial == 0,
                f't1_na rows carry {partial} observed T1/cadd values; policy=fail')
        self.t1_partial_blanked = partial

        want = set(chroms) if chroms else None
        keep = (np.array([c in want for c in chrs], dtype=bool) if want is not None
                else np.ones(len(keys), bool))
        require(keep.any(), 'no cache rows for the requested chromosomes')
        self.X, self.keys, self.genes, self.chrs = X[keep], keys[keep], genes[keep], chrs[keep]
        self.use_chr = sorted(set(self.chrs), key=int)
        require(not chroms or self.use_chr == sorted(set(chroms), key=int),
                'requested chromosome missing from cache')

        self.gl = sorted(set(self.genes))
        gix = {g: i for i, g in enumerate(self.gl)}
        self.gid = np.array([gix[g] for g in self.genes], dtype=np.int32)
        self.ng = len(self.gl)
        gcs = {}
        for g, ch in zip(self.genes, self.chrs):
            gcs.setdefault(g, set()).add(ch)
        require(all(len(v) == 1 for v in gcs.values()), 'multichromosome gene identifiers')
        self.gchr = {g: next(iter(v)) for g, v in gcs.items()}
        order = np.argsort(self.gid, kind='stable')
        cuts = np.searchsorted(self.gid[order], np.arange(self.ng + 1))
        self.pairs_of = [order[cuts[g]:cuts[g + 1]] for g in range(self.ng)]
        require(sum(map(len, self.pairs_of)) == len(self.keys), 'gene pair partition')
        n_genes = self.col('n_genes')
        require(np.isfinite(n_genes).all() and (n_genes >= 1).all(), 'invalid n_genes')
        self.inv_n_genes = (1.0 / n_genes).astype(np.float32)
        _, first = np.unique(self.keys, return_index=True)
        self.unique_rows = np.sort(first)

    def _fp(self, path):
        st = Path(path).stat()
        self.fingerprints[str(path)] = dict(size=st.st_size, mtime_ns=st.st_mtime_ns)

    def col(self, name):
        return self.X[:, self.cols.index(name)].astype(np.float64)

    def pair_weight(self, mode):
        require(mode in ('none', 'inverse-n-genes'), 'unknown --pair-weight')
        return np.ones(len(self.keys), np.float32) if mode == 'none' else self.inv_n_genes

class Dosage:

    def __init__(self, root, chrom, samples=None):
        path = Path(root) / f'annot/ds/chr{chrom}.ds.npz'
        require(path.is_file(), f'chr{chrom} dosage cache missing')
        st = path.stat()
        self.fingerprint = {str(path): dict(size=st.st_size, mtime_ns=st.st_mtime_ns)}
        with np.load(path, allow_pickle=True) as z:
            self.samples = z['samples'].astype(str).tolist()
            require(len(self.samples) == COHORT_N, 'sample count differs from cohort')
            require(len(set(map(norm_iid, self.samples))) == len(self.samples),
                    'normalised sample key duplicates')
            if samples is not None:
                require(self.samples == samples, 'chromosome sample order mismatch')
            indptr, indices, values = z['indptr'], z['indices'], z['data']
            vkeys = z['keys'].astype(str)
            require(values.dtype == np.float32, 'DS data must be prebuilt float32')
            require(len(indices) == len(values) == int(indptr[-1]), 'DS nnz mismatch')
            require(len(indptr) == len(vkeys) + 1, 'DS indptr length')
            require(indptr[0] == 0 and np.all(indptr[1:] >= indptr[:-1]), 'DS indptr monotonicity')
            self.ns = len(self.samples)
            require(indices.min() >= 0 and indices.max() < self.ns, 'DS column bounds')
            self.M = sp.csr_matrix((values, indices.astype(np.int32, copy=False),
                                    indptr.astype(np.int32, copy=False)),
                                   shape=(len(vkeys), self.ns), copy=False)
        self.row_of = {k: i for i, k in enumerate(vkeys)}
        require(len(self.row_of) == len(vkeys), 'DS variant key duplicates')
        self.n_variants = len(vkeys)
        self.nnz = int(self.M.nnz)

    def locate(self, keys):
        missing = int(sum(k not in self.row_of for k in keys))
        require(missing == 0, f'{missing} annotation keys absent from dosage cache')
        return np.array([self.row_of[k] for k in keys], dtype=np.int32)

def burden(dosage, vloc, weight):
    return dosage.M[vloc].T @ weight

def zrow(S):
    Z = np.asarray(S, dtype=np.float64)
    Z = Z - Z.mean(axis=-1, keepdims=True)
    sd = np.sqrt(np.mean(Z * Z, axis=-1, keepdims=True) + NUM_EPS ** 2)
    return Z / sd, sd

def ecdf_avg_rank(values, reference_rows, n_reference):
    ref = np.asarray(values)[reference_rows]
    observed = np.sort(ref[np.isfinite(ref)])
    require(len(observed) > 0, 'annotation has no observed reference values')
    p = np.full(len(values), .5)
    mask = np.isfinite(values)
    left = np.searchsorted(observed, values[mask], side='left')
    right = np.searchsorted(observed, values[mask], side='right')
    p[mask] = np.where(right > left, (left + right + 1) / 2, right) / n_reference
    return p, dict(M=int(n_reference), n_observed=int(len(observed)),
                   n_missing=int((~np.isfinite(ref)).sum()), tie_method='average', missing_weight=.5)

def ecdf_reference(cache, mode):
    require(mode in ('band-all', 'band-chrom', 'band-cadd-available'), 'unknown --ecdf-reference')
    rows = cache.unique_rows
    if mode == 'band-cadd-available':
        have = set(cache.cadd_available)
        rows = rows[np.array([c in have for c in cache.chrs[rows]], dtype=bool)]
    require(len(rows) > 0, 'empty ECDF reference')
    return rows, len(rows)

def phi_arm(cache, arm, v10mod, phi_json, ecdf_mode):
    require(arm in ('learned', 'flat', 'cadd_rank', 'pibar'), 'unknown --arm')
    if arm == 'flat':
        return np.ones(len(cache.keys), np.float64), dict(arm='flat', definition='phi=1')
    if arm == 'learned':
        require(phi_json is not None and v10mod is not None,
                '--phi-json and --v10-script required for the learned arm')
        with open(phi_json) as fh:
            model = json.load(fh)
        require(model.get('version') == 'v10', 'phi export version')
        require(list(model['design']['cols']) == cache.cols,
                'phi design columns differ from the assembled 34-column cache')
        need = {model['design']['rules'][n]['source'] for n in model['design']['order']}
        for cname in sorted(need):
            j = cache.cols.index(cname)
            require(bool(np.isfinite(cache.X[:, j]).any()),
                    f'phi source column {cname} is entirely missing in this band subset')
        v10mod.load_libraries(1)
        phi = np.asarray(v10mod.PhiPredictor(model).predict(cache.X, cache.cols), dtype=np.float64)
        require(np.isfinite(phi).all() and (phi > 0).all() and (phi < 1).all(),
                'phi outside the open unit interval')
        return phi, dict(arm='learned', phi_json_sha256=sha_file(phi_json),
                         intercept=model['default_intercept'],
                         intercept_rule=model['default_intercept_rule'],
                         code_sha256=model.get('code_sha256'),
                         scale='sigmoid PIP score in (0,1)',
                         note='Z_g standardisation makes T invariant to a multiplicative rescaling of '
                              'phi but NOT to an additive shift')
    rows, M = ecdf_reference(cache, ecdf_mode)
    if arm == 'cadd_rank':
        vals = cache.col('cadd')
        require(np.isfinite(vals).any(), 'cadd entirely missing for this subset — extract cadd first')
        p, meta = ecdf_avg_rank(vals, rows, M)
        return p, dict(arm='cadd_rank', definition='pi_hat = average-rank ECDF(cadd)',
                       ecdf=ecdf_mode, **meta)
    parts, metas = [], {}
    for name in STAAR_ANNOTATIONS:
        require(name in cache.cols, f'v9-eight column {name} absent')
        p, m = ecdf_avg_rank(cache.col(name), rows, M)
        parts.append(p)
        metas[name] = m
    pibar = np.mean(np.vstack(parts), axis=0)
    return pibar, dict(arm='pibar',
                       definition='arithmetic mean of the 8 v9-eight average-rank ECDF scores '
                                  '(NOT a Cauchy combination — Cauchy applies to p-values only)',
                       columns=STAAR_ANNOTATIONS, ecdf=ecdf_mode, per_column=metas)

class Residuals:

    def __init__(self, resid_dir, traits, samples, min_n, shuffle_seed, strata_path=None):
        self.names = list(traits)
        root = Path(resid_dir)
        summary_path = root / 'build_resid_v8.summary.json'
        with open(summary_path) as fh:
            self.summary = json.load(fh)
        with open(root / 'build_resid_v8.done') as fh:
            marker = json.load(fh)
        require(marker.get('status') == 'complete', 'residual builder not complete')
        require(marker['summary_sha256'] == sha_file(summary_path), 'residual completion marker')
        require(self.summary['sample_order_sha256'] == digest(samples),
                'builder/genotype sample order mismatch')
        ns = len(samples)
        index = {norm_iid(s): i for i, s in enumerate(samples)}
        self.full = np.full((len(self.names), ns), np.nan, np.float32)
        self.hashes, self.kind, self.matched = {}, {}, {}
        for it, trait in enumerate(self.names):
            path = root / f'{trait}.resid.tsv'
            self.hashes[trait] = sha_file(path)
            require(self.hashes[trait] == self.summary['traits'][trait]['residual_sha256'],
                    f'{trait}: residual content differs from builder manifest')
            self.kind[trait] = self.summary['traits'][trait].get('resid_kind')
            seen, matched = set(), 0
            with open(path, newline='') as fh:
                reader = csv.DictReader(fh, delimiter='\t')
                require({'IID', 'resid'} <= set(reader.fieldnames or []), 'residual schema mismatch')
                for rec in reader:
                    key = norm_iid(rec['IID'])
                    require(key not in seen, 'duplicate normalised residual IID')
                    seen.add(key)
                    if key not in index:
                        continue
                    v = float(rec['resid'])
                    if math.isfinite(v):
                        self.full[it, index[key]] = v
                        matched += 1
            require(matched >= min_n, f'{trait}: insufficient matched residuals')
            self.matched[trait] = matched

        self.shuffled = shuffle_seed is not None
        if self.shuffled:
            perm = np.random.default_rng(shuffle_seed).permutation(ns)
            self.full = self.full[:, perm]
        self.mask = np.isfinite(self.full)
        self.R = np.zeros_like(self.full, dtype=np.float32)
        self.n = self.mask.sum(axis=1)
        for it in range(len(self.names)):
            obs = self.full[it, self.mask[it]].astype(np.float64)
            sd = obs.std(ddof=0)
            require(sd > NUM_EPS, 'zero residual variance')
            self.R[it, self.mask[it]] = (obs - obs.mean()) / sd / np.sqrt(self.n[it])
        strata = [''] * ns
        if strata_path:
            seen = set()
            with open(strata_path, newline='') as fh:
                reader = csv.DictReader(fh, delimiter='\t')
                require({'IID', 'stratum'} <= set(reader.fieldnames or []), 'permutation strata schema')
                for rec in reader:
                    key = norm_iid(rec['IID'])
                    require(key not in seen, 'duplicate permutation-stratum IID')
                    seen.add(key)
                    if key in index:
                        require(bool(rec['stratum']), 'empty stratum')
                        strata[index[key]] = rec['stratum']
            require(all(strata), 'strata missing for genotype samples')
            if self.shuffled:
                strata = [strata[i] for i in perm]
        packed = np.packbits(self.mask.T, axis=1)
        groups = {}
        for i, bits in enumerate(packed):
            groups.setdefault((bits.tobytes(), strata[i]), []).append(i)
        self.groups = [np.asarray(v, np.int32) for v in groups.values()]
        movable = sum(len(v) for v in self.groups if len(v) > 1 and self.mask[:, v[0]].any())
        require(movable >= min_n, 'too few exchangeable observed individuals')
        self.meta = dict(kind='joint within exact observed-trait mask and optional stratum',
                         n_blocks=len(self.groups), movable_observed_samples=int(movable),
                         singleton_samples=int(sum(len(v) == 1 for v in self.groups)),
                         residual_shuffled=self.shuffled, shuffle_seed=shuffle_seed,
                         matched={k: int(v) for k, v in self.matched.items()},
                         n={t: int(x) for t, x in zip(self.names, self.n)}, resid_kind=self.kind,
                         residual_sha256=self.hashes)

    def permutations(self, count, seed):
        rng = np.random.default_rng(seed)
        out = []
        for _ in range(count):
            perm = np.arange(self.R.shape[1], dtype=np.int32)
            for group in self.groups:
                perm[group] = rng.permutation(group)
            out.append(perm)
        return out
