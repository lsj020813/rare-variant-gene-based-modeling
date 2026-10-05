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


import bisect
import collections
import csv
import datetime
import gzip
import io
import json
import math
import os
from pathlib import Path
import random
import re
import resource
import shutil
import hashlib
import fcntl
import threading
import sys
import time
import zipfile

SAMPLE_N = 10000
SEED = 20260928
INTERNAL_KNOTS = 4
KNOT_QUANTILES = (0.2, 0.4, 0.6, 0.8)
MIN_TEAMS = 3
PALINDROMIC_POLICY = "exclude_all"
DS_MISSING_POLICY = "mean_impute_per_variant"
FEATURE_MISSING_POLICY = "observed_member_mean_then_pairwise_complete"
R2_QC = 0.3
LD_R2 = 0.8
LD_DISTANCE_CUT = 0.2
TEAM_BLOCK_SIZE = 2000
FLANK_BP = 3000
VARIANCE_RATIO_MIN = 0.5
PERMUTATIONS = 1000
FDR_ALPHA = 0.10
THREADS = 4
GIB = 1024 ** 3
AS_CAP = 36 * GIB
QUERY_BP = 1000000
DOT_BLOCK = 256
CHUNK_BP = 2000000
BOUNDARY_WINDOW = 200
PROGRESS_GENES = 50
ROOT = Path(_config_path('${PROJECT_ROOT}/work'))
REF = ROOT / 'ref'
GTF = REF / 'deductive/gencode.nochr.gtf.gz'
GPN = REF / 'gpnmsa/scores.tsv.bgz'
BBJ = REF / 'bbj/hum0197.v3.BBJ.Hei.v1.zip'
BBJ_MEMBER = 'hum0197.v3.BBJ.Hei.v1/GWASsummary_Height_Japanese_SakaueKanai2020.auto.txt.gz'
LEADS = ROOT / 'prs/out/hei/lead_snps.tsv'
AUDIT = collections.Counter()
np = None
STOP = threading.Event()
BASE = ROOT / 'phi_gate'
PRIVATE = BASE / 'private'
LOCK_HANDLE = None
CONTRACT = 'phi_gate_v6_r3_stream_block_complete_v1'

class GateError(Exception):
    pass

def need(condition, message):
    if not condition:
        raise GateError(message)

def log(stage, **counts):
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024
    record = dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  stage=stage, rss_peak_gib=round(rss / GIB, 3))
    try:
        pages = int(Path('/proc/self/statm').read_text().split()[1])
        record['rss_current_gib'] = round(pages * os.sysconf('SC_PAGE_SIZE') / GIB, 3)
    except (OSError, ValueError, IndexError):
        record['rss_current_gib'] = None
    record.update(counts)
    print(json.dumps(record, allow_nan=False, sort_keys=True), flush=True)

def chrom(value):
    value = str(value).strip()
    if value.lower().startswith('chr'):
        value = value[3:]
    return str(int(value)) if value.isdigit() else value

def key(value):
    fields = value.strip().split(':')
    need(len(fields) == 4 and fields[1].isdigit(), 'Invalid hg19 allele key.')
    c, p, r, a = chrom(fields[0]), int(fields[1]), fields[2].upper(), fields[3].upper()
    need(p > 0 and r != a and re.fullmatch('[ACGTN]+', r) is not None
         and re.fullmatch('[ACGTN]+', a) is not None,
         'Unsplit multiallelic/symbolic or invalid universe key; allele-specific MAF required.')
    return c, p, r, a

def rc(value):
    return value.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]

def orientation(source_ref, source_alt, target_ref, target_alt):
    for r, a, sign, label in (
            (source_ref, source_alt, 1, 'direct'),
            (source_alt, source_ref, -1, 'swap'),
            (rc(source_ref), rc(source_alt), 1, 'reverse_complement'),
            (rc(source_alt), rc(source_ref), -1, 'swap_reverse_complement')):
        if (r, a) == (target_ref, target_alt):
            return sign, label
    return None

def owned(path):
    need(path.is_dir() and path.stat().st_uid == os.getuid()
         and str(path.resolve()).startswith(_config_path('${PROJECT_ROOT}/')),
         _config_path('Owned directory under ${PROJECT_ROOT} required.'))

def limit_resources():
    sys.dont_write_bytecode = True
    os.umask(0o077)
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    cap = min([AS_CAP] + [x for x in (soft, hard) if x != resource.RLIM_INFINITY])
    resource.setrlimit(resource.RLIMIT_AS, (cap, cap))
    os.sched_setaffinity(0, set(sorted(os.sched_getaffinity(0))[:THREADS]))
    for name in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS',
                 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'ARROW_NUM_THREADS'):
        os.environ[name] = '1'
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    return cap

def process_snapshot():
    records = {}
    heavy = ('deepsea', 'deepripe', 'sparse_genotype', 'sparse-genotype',
             'step2_spa', 'saige_step2', 'saige-step2', 'step2_tests',
             'deeprvat_associate', 'deeprvat-associate', 'annotation', 'annotate')
    for entry in Path('/proc').iterdir():
        if not entry.name.isdigit() or int(entry.name) == os.getpid():
            continue
        try:
            uid = entry.stat().st_uid
            raw = (entry / 'stat').read_text()
            fields = raw[raw.rfind(')') + 2:].split()
            records[(entry.name, fields[19])] = (
                uid, int(fields[11]) + int(fields[12]),
                max(0, int(fields[21])) * os.sysconf('SC_PAGE_SIZE'))
            tokens = [os.path.basename(t).lower()
                      for t in (entry / 'cmdline').read_text().split('\0') if t]
            scripts = [t for t in tokens if t.endswith(('.py', '.r', '.sh'))]
            active = any(h in t for h in heavy for t in scripts)
            active = active or ('deeprvat' in tokens and
                                any(t in tokens for t in ('associate', 'association')))
            need(not active, 'Resource guard: incompatible heavyweight stage active.')
        except (FileNotFoundError, ProcessLookupError):
            continue
    return records

def guard():
    before = process_snapshot()
    start = time.monotonic()
    time.sleep(1)
    after = process_snapshot()
    elapsed = time.monotonic() - start
    cpu, ram = collections.Counter(), collections.Counter()
    for ident, (uid, ticks, rss) in after.items():
        if uid == os.getuid():
            continue
        ram[uid] += rss
        old = before.get(ident, (uid, ticks, 0))[1]
        cpu[uid] += max(0, ticks - old) / os.sysconf('SC_CLK_TCK') / elapsed * 100
    need(max(ram.values(), default=0) <= 32 * GIB and sum(ram.values()) <= 64 * GIB,
         'Resource guard: other-user RAM limit exceeded.')
    need(max(cpu.values(), default=0) <= 400 and sum(cpu.values()) <= 800,
         'Resource guard: other-user CPU limit exceeded.')
    mem = dict(line.split(':', 1) for line in Path('/proc/meminfo').read_text().splitlines())
    need(int(mem['MemAvailable'].split()[0]) * 1024 >= 80 * GIB,
         'Resource guard: available RAM below 80 GiB.')
    need(shutil.disk_usage(ROOT).free >= GIB, 'Resource guard: insufficient output space.')
    log('resource_guard_pass', other_users=len(ram))

def read_universe(ch):
    log('universe_start')
    records = {}
    with gzip.open(ROOT / ('fset/out/uni_chr%s.tsv.gz' % ch), 'rt') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        need(reader.fieldnames is not None and
             {'key', 'pos38', 'maf', 'r2', 'is_coding'} <= set(reader.fieldnames),
             'Universe schema mismatch.')
        for row in reader:
            AUDIT['universe_rows'] += 1
            if int(row['is_coding']) != 0 or float(row['r2']) < R2_QC:
                continue
            need(math.isfinite(float(row['r2'])), 'Invalid universe R2.')
            k = key(row['key'])
            need(k[0] == ch and k not in records, 'Universe chromosome/duplicate allele key error.')
            p, maf = int(row['pos38']), float(row['maf'])
            need(p > 0 and math.isfinite(maf) and 0 <= maf <= 0.5, 'Invalid universe position/MAF.')
            records[k] = (p, maf)
    need(records, 'Empty filtered universe.')
    keys = sorted(records)
    positions = np.array([records[k][0] for k in keys], dtype=np.int64)
    mafs = np.array([records[k][1] for k in keys], dtype=np.float64)
    index = {k: i for i, k in enumerate(keys)}
    AUDIT['universe_filtered'] = len(keys)
    AUDIT['universe_duplicate_keys'] = 0
    log('universe_done', variants=len(keys))
    return keys, positions, mafs, index

def read_leads(ch):
    leads = {}
    with open(LEADS) as handle:
        for line in handle:
            row = line.split()
            need(len(row) == 5, 'Lead cache must be headerless with five columns.')
            if chrom(row[0]) != ch:
                continue
            k = key(row[3])
            need(k[0] == ch and int(row[1]) == k[1] and k not in leads,
                 'Lead hg19 key/position or uniqueness mismatch.')
            need(int(row[4]) > 0, 'Invalid lifted lead position.')
            leads[k] = int(row[4])
    need(leads, 'No same-chromosome leads; distance adjustment undefined.')
    AUDIT['lead_duplicate_keys'] = 0
    return leads

def read_lifts(ch, keys, positions, index, leads):
    log('liftover_start', variants=len(keys), leads=len(leads))
    alleles = [None] * len(keys)
    lead_seen = set()
    with gzip.open(REF / ('lift38_keyed/chr%s.keyed38.vcf.gz' % ch), 'rt') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            row = line.rstrip('\n').split('\t', 8)
            need(len(row) >= 8, 'Invalid keyed38 VCF row.')
            fields = row[2].split(':')
            if len(fields) != 4 or not fields[1].isdigit():
                AUDIT['lift_unparseable_unrelated_ids'] += 1
                continue
            k = (chrom(fields[0]), int(fields[1]), fields[2].upper(), fields[3].upper())
            if k not in index and k not in leads:
                continue
            need(chrom(row[0]) == ch, 'Relevant keyed38 lift moves off chromosome.')
            pos = int(row[1])
            if k in leads:
                need(k not in lead_seen and pos == leads[k], 'Lead lift duplicate/coordinate mismatch.')
                lead_seen.add(k)
            if k in index:
                i = index[k]
                need(alleles[i] is None and pos == int(positions[i]),
                     'Universe keyed38 duplicate/coordinate mismatch.')
                alts = row[4].upper().split(',')
                matches = [(row[3].upper(), a, orientation(k[2], k[3], row[3].upper(), a))
                           for a in alts]
                matches = [m for m in matches if m[2] is not None]
                need(len(matches) == 1, 'Ambiguous or unmatched lifted allele.')
                r, a, orient = matches[0]
                alleles[i] = (r, a, orient[0])
                AUDIT['lift_' + orient[1]] += 1
    need(all(a is not None for a in alleles) and len(lead_seen) == len(leads),
         'Incomplete keyed38 join; no partial analysis permitted.')
    AUDIT['lift_duplicate_keys'] = 0
    log('liftover_done', matched_variants=len(alleles), matched_leads=len(lead_seen))
    return alleles

def assign_genes(ch, positions):
    order = np.argsort(positions, kind='stable')
    ps = positions[order]
    genes = {}
    with gzip.open(GTF, 'rt') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            row = line.rstrip('\n').split('\t')
            if row[2] != 'gene' or chrom(row[0]) != ch:
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]*)"', row[8]))
            if attrs.get('gene_type') != 'protein_coding':
                continue
            gid = attrs['gene_id'].split('.')[0]
            need(gid not in genes, 'Duplicate version-stripped GTF gene ID.')
            lo, hi = int(row[3]), int(row[4])
            need(0 < lo <= hi, 'Invalid GTF interval.')
            genes[gid] = (max(1, lo - FLANK_BP), hi + FLANK_BP)
    assigned = []
    for gid in sorted(genes):
        lo, hi = genes[gid]
        ids = np.sort(order[np.searchsorted(ps, lo, side='left'):
                           np.searchsorted(ps, hi, side='right')])
        if len(ids):
            assigned.append(ids)
    need(assigned, 'No gene assignments.')
    AUDIT['genes_assigned'] = len(assigned)
    AUDIT['gene_variant_instances'] = sum(map(len, assigned))
    return assigned

def load_bed(path, chrom_col, start_col, end_col, val_col, ch, compressed=True):
    intervals = []
    opener = gzip.open if compressed else open
    with opener(path, 'rt') as handle:
        for line in handle:
            row = line.rstrip('\n').split('\t')
            if len(row) <= max(chrom_col, start_col, end_col, val_col or 0):
                continue
            if chrom(row[chrom_col]) != ch:
                continue
            intervals.append((int(row[start_col]), int(row[end_col]),
                              row[val_col] if val_col is not None else '1'))
    intervals.sort()
    need(intervals, 'C-set interval source has no chromosome coverage.')
    return intervals

def overlap_and_dist(intervals, positions):
    starts = [a for a, _, _ in intervals]
    maxend, running = [], -1
    for _, b, _ in intervals:
        running = max(running, b)
        maxend.append(running)
    values, distances = [], np.empty(len(positions), dtype=np.float64)
    for n, p in enumerate(positions):
        p0 = int(p) - 1
        i = bisect.bisect_right(starts, p0) - 1
        j, val = i, None
        while j >= 0 and maxend[j] > p0:
            a, b, v = intervals[j]
            if a <= p0 < b:
                val = v
                break
            j -= 1
        if val is None:
            candidates = []
            if i >= 0:
                candidates.append(p0 - maxend[i] + 1)
            nxt = bisect.bisect_right(starts, p0)
            if nxt < len(intervals):
                candidates.append(intervals[nxt][0] - p0)
            distance = min(candidates)
        else:
            distance = 0
        values.append(val)
        distances[n] = distance
    return values, distances

def annotations(ch, keys, positions, lifted):
    import pyBigWig
    import pysam
    log('annotations_start', variants=len(keys))
    columns = {}
    for name, path, spec, compressed in (
            ('ccre', REF / 'b6_cards/ccre.s.bed', (0, 1, 2, 3), False),
            ('rep', REF / 'repeats/rmsk.txt.gz', (5, 6, 7, 11), True),
            ('cpg', REF / 'cpg/cpgIslandExt.txt.gz', (1, 2, 3, None), True)):
        iv = load_bed(path, *spec, ch=ch, compressed=compressed)
        vals, distances = overlap_and_dist(iv, positions)
        columns[name + '_dist'] = distances
        if name == 'cpg':
            columns['cpg'] = np.array([v is not None for v in vals], dtype=np.float64)
        else:
            vals = np.array([v or 'none' for v in vals], dtype=object)
            for level in sorted(set(vals)):
                columns[name + '_class=' + level] = (vals == level).astype(np.float64)
        del iv, vals
        log('annotations_' + name, variants=len(keys))
    with open(ROOT / 'fset/out/remap_top_tfs.txt') as handle:
        top = [line.strip() for line in handle if line.strip()]
    need(len(top) == 60 and len(set(top)) == 60 and 'n' not in top, 'Expected 60 unique top TF names.')
    perchr = ROOT / ('fset/out/remap_by_chr/chr%s.bed' % ch)
    source = perchr if perchr.exists() else REF / 'remap/remap2022_nr_macs2_hg38_v1_0.bed.gz'
    iv = load_bed(source, 0, 1, 2, 3, ch, compressed=not perchr.exists())
    iv = sorted((a, b, v.split(':')[0]) for a, b, v in iv)
    starts = [v[0] for v in iv]
    maxend = np.maximum.accumulate([v[1] for v in iv])
    for tf in top:
        columns['tf_' + tf] = np.zeros(len(keys), dtype=np.float64)
    columns['tf_n'] = np.zeros(len(keys), dtype=np.float64)
    for n, p in enumerate(positions):
        p0 = int(p) - 1
        j = bisect.bisect_right(starts, p0) - 1
        seen = set()
        while j >= 0 and maxend[j] > p0:
            a, b, tf = iv[j]
            if a <= p0 < b:
                seen.add(tf)
            j -= 1
        columns['tf_n'][n] = len(seen)
        for tf in seen.intersection(top):
            columns['tf_' + tf][n] = 1
    del iv, starts, maxend
    log('annotations_tf', variants=len(keys), tf_columns=len(top))
    mapping = np.full(len(keys), np.nan)
    with pyBigWig.open(str(REF / 'mappability/k36.Umap.MultiTrackMappability.bw')) as bw:
        contigs = [c for c in bw.chroms() if chrom(c) == ch]
        need(len(contigs) == 1, 'BigWig contig absent or ambiguous.')
        contig = contigs[0]
        order = np.argsort(positions, kind='stable')
        blocks = collections.defaultdict(list)
        for i in order:
            blocks[(int(positions[i]) - 1) // QUERY_BP].append(int(i))
        for ids in blocks.values():
            start, end = int(positions[ids[0]]) - 1, int(positions[ids[-1]])
            need(end <= bw.chroms(contig), 'Annotation position outside BigWig chromosome.')
            values = np.asarray(bw.values(contig, start, end, numpy=True))
            mapping[ids] = values[positions[ids] - 1 - start]
    need(np.isfinite(mapping).any(), 'Zero mappability coverage.')
    columns['map_k36'] = mapping
    log('annotations_mappability', matched=int(np.isfinite(mapping).sum()))
    scores = np.full(len(keys), np.nan)
    byblock = collections.defaultdict(list)
    for i, p in enumerate(positions):
        byblock[(int(p) - 1) // QUERY_BP].append(i)
    with pysam.TabixFile(str(GPN), threads=1) as tab:
        contigs = [c for c in tab.contigs if chrom(c) == ch]
        need(len(contigs) == 1, 'GPN tabix contig absent or ambiguous.')
        for block_no, ids in enumerate(byblock.values(), 1):
            bypos = collections.defaultdict(list)
            for i in ids:
                bypos[int(positions[i])].append(i)
            records, seen_gpn = {}, set()
            for line in tab.fetch(contigs[0], min(bypos) - 1, max(bypos)):
                row = line.split('\t')
                need(len(row) == 5 and chrom(row[0]) == ch, 'GPN five-column schema mismatch.')
                p = int(row[1])
                if p not in bypos:
                    continue
                k = (p, row[2].upper(), row[3].upper())
                need(k not in seen_gpn, 'Duplicate GPN allele key.')
                seen_gpn.add(k)
                value = float(row[4])
                if math.isfinite(value):
                    records[k] = value
            for p, local_ids in bypos.items():
                for i in local_ids:
                    r, a, lift_sign = lifted[i]
                    if len(r) != 1 or len(a) != 1 or r not in 'ACGT' or a not in 'ACGT':
                        AUDIT['gpn_indel_or_non_snv'] += 1
                        continue
                    for rr, aa, sign in ((r, a, 1), (rc(r), rc(a), 1),
                                         (a, r, -1), (rc(a), rc(r), -1)):
                        found = records.get((p, rr, aa))
                        if found is not None:
                            scores[i] = found * sign * lift_sign
                            AUDIT['gpn_matched'] += 1
                            break
            if block_no % 5 == 0:
                log('annotations_gpn_progress', blocks_done=block_no, blocks_total=len(byblock))
    need(np.isfinite(scores).any(), 'Zero allele-aware GPN join.')
    columns['gpn_score'] = scores
    names = sorted(columns)
    need(not any(n in names for n in ('maf', 'r2', 'typed', 'team_size')), 'Forbidden feature.')
    matrix = np.column_stack([columns[n] for n in names])
    coverage = {n: dict(observed=int(np.isfinite(columns[n]).sum()), total=len(keys)) for n in names}
    log('annotations_done', features=len(names), variants=len(keys))
    return names, matrix, coverage

def match_bbj(ch, keys, index):
    log('bbj_start', variants=len(keys))
    beta = np.full(len(keys), np.nan)
    matched = np.zeros(len(keys), dtype=bool)
    with zipfile.ZipFile(BBJ) as archive:
        with archive.open(BBJ_MEMBER) as member:
            with gzip.GzipFile(fileobj=member) as zipped:
                with io.TextIOWrapper(zipped) as handle:
                    header = handle.readline().split()
                    need(len(header) == len(set(header)) and
                         {'CHR', 'BP', 'ALLELE1', 'ALLELE0', 'BETA'} <= set(header), 'BBJ schema mismatch.')
                    ix = {n: header.index(n) for n in ('CHR', 'BP', 'ALLELE1', 'ALLELE0', 'BETA')}
                    for line_no, line in enumerate(handle, 1):
                        row = line.split()
                        need(len(row) == len(header), 'Invalid BBJ row length.')
                        if line_no % 500000 == 0:
                            log('bbj_progress', source_rows=line_no, chromosome_rows=AUDIT['bbj_chr_rows'],
                                matched=int(matched.sum()))
                        if chrom(row[ix['CHR']]) != ch:
                            continue
                        AUDIT['bbj_chr_rows'] += 1
                        p = int(row[ix['BP']])
                        a0, a1 = row[ix['ALLELE0']].upper(), row[ix['ALLELE1']].upper()
                        pal = len(a0) == len(a1) == 1 and {a0, a1} in ({'A', 'T'}, {'C', 'G'})
                        if pal:
                            AUDIT['bbj_palindromic_rows_excluded'] += 1
                            continue
                        alternatives = [(a0, a1, 1, 'direct'), (a1, a0, -1, 'swap')]
                        if len(a0) == len(a1) == 1 and a0 in 'ACGT' and a1 in 'ACGT':
                            alternatives += [(rc(a0), rc(a1), 1, 'reverse_complement'),
                                             (rc(a1), rc(a0), -1, 'swap_reverse_complement')]
                        hits = {}
                        for r, a, sign, label in alternatives:
                            i = index.get((ch, p, r, a))
                            if i is not None:
                                hits[i] = (sign, label)
                        need(len(hits) <= 1, 'Ambiguous BBJ match to multiple universe allele keys.')
                        if not hits:
                            continue
                        i, (sign, label) = next(iter(hits.items()))
                        need(not matched[i], 'Duplicate BBJ allele mapping; no arbitrary row selection.')
                        matched[i] = True
                        AUDIT['bbj_' + label] += 1
                        try:
                            effect = float(row[ix['BETA']])
                        except ValueError:
                            AUDIT['bbj_invalid_beta'] += 1
                            continue
                        if not math.isfinite(effect):
                            AUDIT['bbj_invalid_beta'] += 1
                            continue
                        aligned = sign * effect
                        beta[i] = abs(aligned)
    need(matched.any(), 'Zero BBJ allele join.')
    AUDIT['bbj_matched_variants'] = int(matched.sum())
    AUDIT['bbj_finite_beta_variants'] = int(np.isfinite(beta).sum())
    AUDIT['bbj_duplicate_keys'] = 0
    log('bbj_done', matched=int(matched.sum()), total=len(keys), matching_rate=float(matched.mean()))
    return beta, matched

def natural_spline(x):
    low, high = float(x.min()), float(x.max())
    need(high > low, 'Spline covariate has no variation.')
    z = (x - low) / (high - low)
    interior = np.quantile(z, KNOT_QUANTILES)
    knots = np.r_[0.0, interior, 1.0]
    need(len(interior) == INTERNAL_KNOTS and np.all(np.diff(knots) > 0),
         'Spline quantile knots coincide; four distinct internal knots required.')
    last, penultimate = knots[-1], knots[-2]
    last_cube = np.maximum(z - last, 0.0) ** 3

    def d(knot):
        return (np.maximum(z - knot, 0.0) ** 3 - last_cube) / (last - knot)

    reference = d(penultimate)
    basis = np.column_stack([z] + [d(k) - reference for k in knots[:-2]])
    return basis

def residual_target(beta, mafs, distance):
    valid = (np.isfinite(beta) & (beta > 0) & np.isfinite(mafs) &
             (mafs > 0) & np.isfinite(distance) & (distance >= 0))
    need(int(valid.sum()) > 11, 'Insufficient valid BBJ targets for spline OLS.')
    x = np.column_stack([np.ones(int(valid.sum())), natural_spline(np.log10(mafs[valid])),
                         natural_spline(np.log10(distance[valid] + 1.0))])
    y = np.log(beta[valid])
    coef, _, rank, singular = np.linalg.lstsq(x, y, rcond=None)
    need(rank == x.shape[1], 'Spline design is rank deficient.')
    residual = y - x.dot(coef)
    need(np.isfinite(residual).all() and np.var(residual) > 0, 'Invalid or constant residual.')
    result = np.full(len(beta), np.nan)
    result[valid] = residual
    log('pooled_residuals_done', fitted=int(valid.sum()), design_rank=int(rank))
    return result, dict(n_fitted=int(valid.sum()), design_rank=int(rank),
                        design_columns=int(x.shape[1]), response='log_abs_beta',
                        condition_number=float(singular[0] / singular[-1]),
                        residual_variance=float(np.var(residual)))

def dosage_header(ch):
    path = REF / ('orig_index/chr%s.vcf.gz' % ch)
    contigs, ds_number = [], None
    samples = None
    with gzip.open(path, 'rt') as handle:
        for line in handle:
            if line.startswith('##contig=<'):
                match = re.search(r'ID=([^,>]+)', line)
                if match and chrom(match.group(1)) == ch:
                    contigs.append(match.group(1))
            if line.startswith('##FORMAT=<ID=DS,'):
                match = re.search(r'Number=([^,>]+)', line)
                ds_number = match.group(1) if match else None
            if line.startswith('#CHROM\t'):
                samples = line.rstrip('\n').split('\t')[9:]
                break
    need(len(contigs) == 1 and samples is not None and ds_number in ('1', 'A'),
         'Dosage header contig/DS schema invalid.')
    need(len(samples) == len(set(samples)) and len(samples) >= SAMPLE_N,
         'Duplicate or fewer than requested VCF samples.')
    need(all(s and '\n' not in s and '\r' not in s for s in samples), 'Invalid sample header.')
    samples = sorted(samples)
    indices = sorted(random.Random(SEED).sample(range(len(samples)), SAMPLE_N))
    chosen = [samples[i] for i in indices]
    log('dosage_header_done', available_samples=len(samples), selected_samples=len(chosen),
        normalized_contig_matches=len(contigs))
    return path, contigs[0], chosen, ds_number

def complete_teams(matrix):
    from scipy.cluster.hierarchy import fcluster, linkage
    n = len(matrix)
    if n == 1:
        return [np.array([0], dtype=np.int64)]
    dist = np.empty(n * (n - 1) // 2, dtype=np.float64)
    for start in range(0, n, DOT_BLOCK):
        stop = min(n, start + DOT_BLOCK)
        corr = matrix[start:stop].dot(matrix.T)
        np.square(corr, out=corr)
        np.clip(corr, 0, 1, out=corr)
        for j in range(start, stop):
            r2 = corr[j - start, j + 1:]
            d = np.where(r2 > LD_R2, 1.0 - r2, np.maximum(1.0 - r2, LD_DISTANCE_CUT))
            offset = n * j - j * (j + 1) // 2
            dist[offset:offset + len(d)] = d
    tree = linkage(dist, method='complete')
    labels = fcluster(tree, t=np.nextafter(LD_DISTANCE_CUT, -np.inf), criterion='distance')
    groups = collections.defaultdict(list)
    for i, label in enumerate(labels):
        groups[int(label)].append(i)
    teams = [np.array(g, dtype=np.int64) for g in sorted(groups.values(), key=lambda g: g[0])]
    for members in teams:
        if len(members) == 1:
            continue
        for lo in range(0, len(members), DOT_BLOCK):
            left = matrix[members[lo:lo + DOT_BLOCK]]
            for hi in range(lo, len(members), DOT_BLOCK):
                right = left if hi == lo else matrix[members[hi:hi + DOT_BLOCK]]
                r2 = np.clip(left.dot(right.T) ** 2, 0, 1)
                if hi == lo:
                    mask = np.triu(np.ones(r2.shape, dtype=bool), 1)
                    need(np.all(r2[mask] > LD_R2), 'Block team clique invariant failed.')
                else:
                    need(np.all(r2 > LD_R2), 'Block team clique invariant failed.')
    need(np.array_equal(np.sort(np.concatenate(teams)), np.arange(n)),
         'Block teams do not partition the block.')
    return teams

def attach_team_targets(grouped, memberships, target):
    for (_, y), teams in zip(grouped, memberships):
        for t, members in enumerate(teams):
            values = target[members]
            valid = np.isfinite(values)
            if valid.any():
                y[t] = np.max(values[valid])
    AUDIT['teams_with_target'] = sum(int(np.isfinite(y).sum()) for _, y in grouped)
    log('team_targets_done', teams_with_target=AUDIT['teams_with_target'], teams=AUDIT['teams'])

def measure_one(grouped, names):
    log('measurement_1_start', features=len(names))
    result = []
    for j, name in enumerate(names):
        groups = [x[np.isfinite(x[:, j]), j] for x, _ in grouped]
        groups = [x for x in groups if len(x)]
        n = sum(len(x) for x in groups)
        if n:
            overall = sum(float(x.sum()) for x in groups) / n
            within = sum(float(np.sum((x - x.mean()) ** 2)) for x in groups)
            total = sum(float(np.sum((x - overall) ** 2)) for x in groups)
        else:
            within, total = 0.0, 0.0
        ratio = within / total if total > 0 else None
        if ratio is not None:
            need(-1e-10 <= ratio <= 1 + 1e-10, 'Variance decomposition failed.')
            ratio = max(0.0, min(1.0, ratio))
        result.append(dict(annotation=name, observed_teams=n, genes_with_observations=len(groups),
                           within_variance=within / n if n else None,
                           total_variance=total / n if n else None,
                           variance_ratio=ratio, candidate=ratio is not None and ratio >= VARIANCE_RATIO_MIN))
    log('measurement_1_done', candidates=sum(r['candidate'] for r in result), features=len(names))
    return result

def bh(pvalues):
    p = np.asarray(pvalues, dtype=float)
    if not len(p):
        return p
    order = np.argsort(p, kind='stable')
    adjusted = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    q = np.empty(len(p))
    q[order] = np.minimum(adjusted, 1.0)
    return q

def measure_two(grouped, names, first):
    from scipy.stats import rankdata
    candidates = [j for j, r in enumerate(first) if r['candidate']]
    log('measurement_2_start', candidates=len(candidates), permutations=PERMUTATIONS)
    results = []
    for feature_no, j in enumerate(candidates, 1):
        rng = np.random.default_rng(np.random.SeedSequence([SEED, j]))
        observed_sum, null_sum, used = 0.0, np.zeros(PERMUTATIONS), 0
        too_few, constant = 0, 0
        for gene_no, (x, y) in enumerate(grouped, 1):
            if gene_no % PROGRESS_GENES == 0:
                log('measurement_2_gene_progress', annotation_number=feature_no,
                    annotations_total=len(candidates), genes_done=gene_no, genes_total=len(grouped))
            valid = np.isfinite(x[:, j]) & np.isfinite(y)
            if int(valid.sum()) < MIN_TEAMS:
                too_few += 1
                continue
            xr, yr = rankdata(x[valid, j], method='average'), rankdata(y[valid], method='average')
            xr -= xr.mean()
            yr -= yr.mean()
            nx, ny = np.linalg.norm(xr), np.linalg.norm(yr)
            if nx == 0 or ny == 0:
                constant += 1
                continue
            xr, yr = xr / nx, yr / ny
            observed_sum += float(xr.dot(yr))
            for b in range(0, PERMUTATIONS, 16):
                perms = np.stack([rng.permutation(yr)
                                  for _ in range(min(16, PERMUTATIONS - b))])
                null_sum[b:b + len(perms)] += perms.dot(xr)
            used += 1
        if used:
            observed = observed_sum / used
            null = null_sum / used
            null95 = float(np.quantile(null, 0.95))
            pvalue = float((1 + np.sum(null >= observed)) / (PERMUTATIONS + 1))
            null_mean = float(null.mean())
        else:
            observed, null95, null_mean, pvalue = None, None, None, 1.0
        results.append(dict(annotation=names[j], eligible_genes=used, too_few_genes=too_few,
                            constant_genes=constant, observed=observed, null_mean=null_mean,
                            null_95=null95, p=pvalue, status='tested' if used else 'not_estimable'))
        log('measurement_2_progress', annotations_done=feature_no, annotations_total=len(candidates),
            eligible_genes=used)
    qvalues = bh([r['p'] for r in results])
    for r, q in zip(results, qvalues):
        r['q'] = float(q)
        r['passed'] = bool(r['observed'] is not None and r['observed'] > r['null_95'] and q < FDR_ALPHA)
    log('measurement_2_done', tested=sum(r['status'] == 'tested' for r in results),
        passed=sum(r['passed'] for r in results))
    return results

def atomic_json(path, payload):
    tmp = path.with_name(path.name + '.tmp')
    need(not tmp.exists() and not path.exists(), 'Output already exists; manual review required.')
    with open(tmp, 'x') as handle:
        json.dump(payload, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write('\n')
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp, path)
    fd = os.open(str(path.parent), os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)

def sample_digest(samples):
    h = hashlib.sha256()
    for sample in samples:
        h.update(sample.encode('utf-8'))
        h.update(b'\0')
    return h.hexdigest()

class StreamingBlocks:
    def __init__(self, n):
        self.n, self.start = n, 0
        self.stop = min(TEAM_BLOCK_SIZE, n)
        self.buffer = np.empty((min(TEAM_BLOCK_SIZE, n), SAMPLE_N), dtype=np.float64)
        self.filled = 0
        self.pending = {}
        self.pending_peak = 0
        self.labels = np.full(n, -1, dtype=np.int64)
        self.team_count, self.blocks_done = 0, 0
        self.previous_tail = None
        self.touched = np.zeros(n, dtype=bool)
        self.boundaries = []
        self.block_sizes = []

    def push(self, i, row):
        need(self.start <= i < min(self.stop + TEAM_BLOCK_SIZE, self.n),
             'Dosage missing/order error or same-position lookahead exceeds one next block.')
        if i < self.stop:
            self.buffer[i - self.start] = row
            self.filled += 1
        else:
            self.pending[i] = row
            self.pending_peak = max(self.pending_peak, len(self.pending))
        while self.start < self.n and self.filled == self.stop - self.start:
            self.finish_block()

    def finish_block(self):
        matrix = self.buffer[:self.stop - self.start]
        teams = complete_teams(matrix)
        for team in teams:
            self.labels[self.start + team] = self.team_count
            self.team_count += 1
        self.block_sizes.append(len(matrix))
        self.blocks_done += 1
        if self.previous_tail is not None:
            head = matrix[:BOUNDARY_WINDOW]
            r2 = np.clip(self.previous_tail.dot(head.T) ** 2, 0, 1)
            hits = r2 > LD_R2
            left, right = hits.any(axis=1), hits.any(axis=0)
            self.touched[self.start - len(left):self.start] |= left
            self.touched[self.start:self.start + len(right)] |= right
            self.boundaries.append(dict(
                boundary_number=len(self.boundaries) + 1,
                left_window_variants=len(left), right_window_variants=len(right),
                tested_pairs=int(hits.size), high_ld_pairs=int(hits.sum()),
                left_involved_variants=int(left.sum()),
                right_involved_variants=int(right.sum()),
                involved_variants=int(left.sum() + right.sum())))
        self.previous_tail = (matrix[-BOUNDARY_WINDOW:].copy()
                              if self.stop < self.n else None)
        self.start = self.stop
        self.stop = min(self.start + TEAM_BLOCK_SIZE, self.n)
        self.filled = 0
        for i in sorted(list(self.pending)):
            if i < self.stop:
                self.buffer[i - self.start] = self.pending.pop(i)
                self.filled += 1
        if self.blocks_done % 5 == 0 or self.start == self.n:
            log('block_teams_progress', blocks_done=self.blocks_done,
                variants_clustered=self.start, variants_total=self.n,
                block_teams=self.team_count)
        if self.start == self.n:
            self.buffer = None

    def finish(self):
        need(self.start == self.n and not self.pending and np.all(self.labels >= 0),
             'Incomplete dosage/block join; no partial teams allowed.')
        need(len(self.boundaries) == max(0, self.blocks_done - 1),
             'Missing adjacent-block boundary diagnostic.')
        return self.labels, dict(
            definition='adjacent_block_last_200_by_first_200',
            diagnostic_only=True, exhaustive_cross_block_scan=False,
            threshold='r_squared_strictly_greater_than_0.8',
            window_variants=BOUNDARY_WINDOW, blocks=self.blocks_done,
            block_sizes=self.block_sizes, boundaries=len(self.boundaries),
            tested_pairs=sum(b['tested_pairs'] for b in self.boundaries),
            high_ld_pairs=sum(b['high_ld_pairs'] for b in self.boundaries),
            involved_unique_variants=int(self.touched.sum()),
            per_boundary=self.boundaries,
            lookahead_peak_variants=self.pending_peak,
            limitations=['Pairs beyond either 200-variant window are not inspected.',
                         'Nonadjacent blocks are not inspected.',
                         'Counts do not establish permutation calibration.'])

def stream_block_teams(ch, keys, index):
    from cyvcf2 import VCF
    need(keys and all(k[0] == ch for k in keys) and
         all(keys[i - 1] < keys[i] for i in range(1, len(keys))),
         'Assigned keys must be unique in hg19 position/full-key order.')
    path, contig, samples, number = dosage_header(ch)
    reader = VCF(str(path), samples=samples, lazy=True, threads=1)
    state = StreamingBlocks(len(keys))
    try:
        actual = reader.samples
        need(len(actual) == SAMPLE_N and set(actual) == set(samples),
             'cyvcf2 sample selection mismatch.')
        lookup = {s: j for j, s in enumerate(actual)}
        reorder = np.array([lookup[s] for s in samples], dtype=np.int64)
        found = np.zeros(len(keys), dtype=bool)
        regions = sorted({(k[1] - 1) // CHUNK_BP for k in keys})
        previous_pos = 0
        for done, region in enumerate(regions, 1):
            lo, hi = region * CHUNK_BP + 1, (region + 1) * CHUNK_BP
            for record in reader('%s:%d-%d' % (contig, lo, hi)):
                if not lo <= record.POS <= hi:
                    continue
                need(chrom(record.CHROM) == ch and record.POS >= previous_pos,
                     'VCF chromosome/position stream order mismatch.')
                previous_pos = record.POS
                alts = [a.upper() for a in record.ALT]
                wanted = [(ai, index.get((ch, record.POS, record.REF.upper(), alt)))
                          for ai, alt in enumerate(alts)]
                wanted = sorted(((ai, i) for ai, i in wanted if i is not None),
                                key=lambda item: item[1])
                if not wanted:
                    continue
                need(len(alts) == 1 or number == 'A', 'Multiallelic DS requires Number=A.')
                ds = record.format('DS')
                need(ds is not None and ds.ndim == 2 and ds.shape == (SAMPLE_N, len(alts)),
                     'Dosage DS dimension/allele count mismatch.')
                for ai, i in wanted:
                    need(not found[i], 'Duplicate dosage allele key.')
                    values = np.array(ds[reorder, ai], dtype=np.float64, copy=True)
                    good = np.isfinite(values)
                    need(good.any() and np.all((values[good] >= 0) & (values[good] <= 2)),
                         'All-missing or invalid dosage.')
                    AUDIT['ds_missing_cells_unique'] += int((~good).sum())
                    values[~good] = values[good].mean()
                    values -= values.mean()
                    norm = float(np.linalg.norm(values))
                    if norm > 0:
                        values /= norm
                    else:
                        values[:] = 0
                        AUDIT['ds_constant_unique_variants'] += 1
                    found[i] = True
                    state.push(i, values)
                del ds, values
            if done % 5 == 0 or done == len(regions):
                log('dosage_progress', regions_done=done, regions_total=len(regions),
                    variants_found=int(found.sum()))
        need(found.all(), 'Incomplete dosage join; no partial teams allowed.')
        AUDIT['dosage_duplicate_keys'] = 0
        labels, boundary = state.finish()
        AUDIT['chromosome_blocks'] = state.blocks_done
        AUDIT['chromosome_block_teams'] = state.team_count
    finally:
        reader.close()
    return labels, sample_digest(samples), boundary

def build_teams(assignments, features, labels):
    grouped, memberships = [], []
    sizes = collections.Counter()
    for g, ids in enumerate(assignments, 1):
        local = collections.defaultdict(list)
        for i in ids:
            local[int(labels[i])].append(int(i))
        members = [np.asarray(local[label], dtype=np.int64) for label in sorted(local)]
        need(sum(map(len, members)) == len(ids) and
             np.array_equal(np.sort(np.concatenate(members)), ids),
             'Teams do not partition the gene.')
        x = np.full((len(members), features.shape[1]), np.nan)
        for t, variants in enumerate(members):
            need(len(variants) <= TEAM_BLOCK_SIZE and
                 np.all(labels[variants] == labels[variants[0]]) and
                 np.all(variants // TEAM_BLOCK_SIZE == variants[0] // TEAM_BLOCK_SIZE),
                 'Gene-team must be a subset of one verified block team.')
            sub = features[variants]
            finite = np.isfinite(sub)
            count = finite.sum(axis=0)
            x[t] = np.divide(np.where(finite, sub, 0).sum(axis=0), count,
                             out=np.full(features.shape[1], np.nan), where=count > 0)
            sizes[len(variants)] += 1
        grouped.append((x, np.full(len(members), np.nan)))
        memberships.append(members)
        if g % PROGRESS_GENES == 0 or g == len(assignments):
            log('teams_progress', genes_done=g, genes_total=len(assignments),
                team_instances=sum(sizes.values()))
    AUDIT['teams'] = sum(sizes.values())
    return grouped, memberships, dict(team_instances=sum(sizes.values()),
                                      singleton_instances=sizes[1], largest_team=max(sizes),
                                      definition='chromosome_block_team_intersect_gene',
                                      block_variants=TEAM_BLOCK_SIZE)

def make_owned_directory(path):
    if path == ROOT:
        owned(path)
        return
    need(str(path).startswith(str(ROOT) + '/'), 'Write outside owned workspace.')
    make_owned_directory(path.parent)
    need(not path.is_symlink(), 'Symlink output path is forbidden.')
    path.mkdir(mode=0o700, exist_ok=True)
    owned(path)
    if path == PRIVATE or PRIVATE in path.parents:
        os.chmod(path, 0o700)

def heartbeat():
    while not STOP.wait(30):
        try:
            guard()
            log('heartbeat')
        except BaseException:
            log('failed', reason='Periodic resource guard failed; no completion marker.')
            os._exit(1)

def digest_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for data in iter(lambda: handle.read(8 * 1024 ** 2), b''):
            h.update(data)
    return h.hexdigest()

def versions():
    from importlib.metadata import version
    return dict(python=sys.version.split()[0],
                **{name: version(name) for name in ('numpy', 'scipy', 'cyvcf2', 'pysam', 'pyBigWig')})

def implementation():
    return dict(contract=CONTRACT, source_sha256=digest_file(Path(__file__)),
                seal_sha256=digest_file(BASE / 'height_seed_phi_prereg_20260928.md'),
                versions=versions(),
                teams=dict(definition='chromosome_block_team_intersect_gene',
                           block_variants=TEAM_BLOCK_SIZE,
                           order='hg19_position_then_normalized_full_key',
                           linkage='scipy_complete', normalized_dtype='float64',
                           cut='nextafter(0.2,-inf)', strict_r_squared=LD_R2,
                           boundary_window=BOUNDARY_WINDOW,
                           diagnostics_affect_decision=False),
                defaults=dict(samples=SAMPLE_N, seed=SEED, min_teams=MIN_TEAMS,
                              knot_quantiles=list(KNOT_QUANTILES), response='log_abs_beta',
                              ld_missing=DS_MISSING_POLICY, ld_imputation_dtype='float64',
                              feature_missing=FEATURE_MISSING_POLICY,
                              palindromic=PALINDROMIC_POLICY))

def input_signature(ch):
    paths = [ROOT / ('fset/out/uni_chr%s.tsv.gz' % ch), GTF, GPN, BBJ, LEADS,
             REF / ('orig_index/chr%s.vcf.gz' % ch),
             REF / ('lift38_keyed/chr%s.keyed38.vcf.gz' % ch),
             REF / 'b6_cards/ccre.s.bed', REF / 'repeats/rmsk.txt.gz',
             REF / 'cpg/cpgIslandExt.txt.gz',
             REF / 'mappability/k36.Umap.MultiTrackMappability.bw',
             ROOT / 'fset/out/remap_top_tfs.txt']
    remap = ROOT / ('fset/out/remap_by_chr/chr%s.bed' % ch)
    paths.append(remap if remap.exists() else REF / 'remap/remap2022_nr_macs2_hg38_v1_0.bed.gz')
    for indexed in (GPN, REF / ('orig_index/chr%s.vcf.gz' % ch)):
        found = [Path(str(indexed) + ext) for ext in ('.tbi', '.csi')
                 if Path(str(indexed) + ext).is_file()]
        need(found, 'Required tabix/VCF index absent.')
        paths.extend(found)
    result = {}
    for path in paths:
        stat = path.stat()
        result[str(path.relative_to(ROOT))] = [stat.st_size, stat.st_mtime_ns,
                                               stat.st_ino, stat.st_dev]
    return result

def fresh_directory(path):
    make_owned_directory(path)
    need(not any(path.iterdir()), 'Existing done or partial output; manual review required.')
    lock = path / '.running'
    lock.mkdir(mode=0o700)
    return lock

def atomic_npz(path, arrays):
    tmp = path.with_name(path.name + '.tmp')
    need(not tmp.exists() and not path.exists(), 'Private artifact already exists.')
    with open(tmp, 'xb') as handle:
        np.savez(handle, **arrays)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp, path)
    fd = os.open(str(path.parent), os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)

def nearest_lead(positions, leads):
    lp = np.array(sorted(set(leads.values())), dtype=np.int64)
    j = np.searchsorted(lp, positions)
    return np.minimum(np.abs(positions - lp[np.maximum(0, j - 1)]),
                      np.abs(positions - lp[np.minimum(len(lp) - 1, j)])).astype(float)

def run_chr(ch):
    out, private = BASE / ('out/chr%s' % ch), PRIVATE / ('chr%s' % ch)
    outlock = fresh_directory(out)
    privlock = fresh_directory(private)
    provenance, signatures = implementation(), input_signature(ch)
    keys, positions, mafs, index = read_universe(ch)
    assignments = assign_genes(ch, positions)
    selected = np.unique(np.concatenate(assignments))
    remap = np.full(len(keys), -1, dtype=np.int64)
    remap[selected] = np.arange(len(selected))
    assignments = [remap[ids] for ids in assignments]
    keys = [keys[int(i)] for i in selected]
    positions, mafs = positions[selected], mafs[selected]
    index = {k: i for i, k in enumerate(keys)}
    del remap, selected
    AUDIT['assigned_unique_variants'] = len(keys)
    leads = read_leads(ch)
    lifted = read_lifts(ch, keys, positions, index, leads)
    names, features, coverage = annotations(ch, keys, positions, lifted)
    del lifted
    largest = max(map(len, assignments))
    artifact_estimate = (sum(map(len, assignments)) * (len(names) * 8 + 32) +
                         len(keys) * 96)
    disk_need = artifact_estimate + 8 * GIB
    need(shutil.disk_usage(private).free >= disk_need,
         'Insufficient private output space; no variant or gene is excluded.')
    log('disk_plan', required_gib=round(disk_need / GIB, 2), largest_gene_variants=largest,
        assigned_variants=len(keys), genes=len(assignments))
    labels, samples_sha, boundary = stream_block_teams(ch, keys, index)
    grouped, memberships, team_summary = build_teams(assignments, features, labels)
    del features, assignments, labels
    descriptive = measure_one(grouped, names)
    for row in descriptive:
        row.pop('candidate')
    beta, matched = match_bbj(ch, keys, index)
    distance = nearest_lead(positions, leads)
    variant_hash = np.empty((len(keys), 32), dtype=np.uint8)
    for i, k in enumerate(keys):
        value = ':'.join(map(str, k)).encode('ascii')
        variant_hash[i] = np.frombuffer(hashlib.sha256(value).digest(), dtype=np.uint8)
    del keys, index, positions, leads
    gene_lengths = np.array([len(x) for x, _ in grouped], dtype=np.int64)
    member_lengths = np.array([len(m) for gene in memberships for m in gene], dtype=np.int64)
    arrays = dict(X=np.vstack([x for x, _ in grouped]), gene_lengths=gene_lengths,
                  member_lengths=member_lengths,
                  member_flat=np.concatenate([m for gene in memberships for m in gene]),
                  maf=mafs, absbeta=beta, lead=distance, variant_sha256=variant_hash,
                  names=np.asarray(names))
    artifact = private / 'chr_build.npz'
    atomic_npz(artifact, arrays)
    need(input_signature(ch) == signatures and implementation() == provenance,
         'Inputs or implementation changed during chromosome build.')
    diagnostics = dict(chromosome=int(ch), descriptive_only=True, provenance=provenance,
                       audit=dict(AUDIT), annotations=coverage, teams=team_summary,
                       boundary_diagnostics=boundary,
                       measurement_1_descriptive=descriptive,
                       bbj_matching_rate=float(matched.mean()),
                       bbj_finite_beta_variants=int(np.isfinite(beta).sum()),
                       bbj_zero_beta_variants=int(np.sum(beta == 0)),
                       variants=len(mafs), genes=len(gene_lengths),
                       sample_sha256=samples_sha,
                       projected_disk_gib=round(disk_need / GIB, 2))
    atomic_json(out / 'diagnostics.json', diagnostics)
    manifest = dict(chromosome=int(ch), provenance=provenance, inputs=signatures,
                    sample_sha256=samples_sha, artifact_sha256=digest_file(artifact),
                    diagnostics_sha256=digest_file(out / 'diagnostics.json'))
    atomic_json(private / 'manifest.json', manifest)
    privlock.rmdir()
    outlock.rmdir()
    atomic_json(out / 'chr.done', dict(status='complete', chromosome=int(ch),
                manifest_sha256=digest_file(private / 'manifest.json'),
                outputs_relative_to_phi_gate=['out/chr%s/diagnostics.json' % ch,
                         'private/chr%s/chr_build.npz' % ch,
                         'private/chr%s/manifest.json' % ch]))
    log('chr_complete', chromosome=int(ch))

def read_built(ch, provenance):
    out, private = BASE / ('out/chr%s' % ch), PRIVATE / ('chr%s' % ch)
    owned(private)
    need(private.stat().st_mode & 0o077 == 0, 'Private directory permissions are not 0700.')
    need(not (out / '.running').exists() and not (private / '.running').exists(),
         'Chromosome build is incomplete.')
    with open(out / 'chr.done') as handle:
        done = json.load(handle)
    with open(private / 'manifest.json') as handle:
        manifest = json.load(handle)
    need(done['status'] == 'complete' and done['chromosome'] == int(ch) and
         done['manifest_sha256'] == digest_file(private / 'manifest.json'),
         'Invalid chromosome completion marker.')
    need(manifest['chromosome'] == int(ch) and manifest['provenance'] == provenance,
         'Mixed implementation/seal/library versions.')
    need(manifest['inputs'] == input_signature(ch), 'Chromosome inputs changed after build.')
    need(manifest['diagnostics_sha256'] == digest_file(out / 'diagnostics.json') and
         manifest['artifact_sha256'] == digest_file(private / 'chr_build.npz'),
         'Completed chromosome artifact checksum mismatch.')
    with open(out / 'diagnostics.json') as handle:
        diagnostics = json.load(handle)
    with np.load(private / 'chr_build.npz', allow_pickle=False) as archive:
        data = {k: archive[k] for k in archive.files}
    n = len(data['maf'])
    need(len(data['absbeta']) == n and len(data['lead']) == n and
         data['variant_sha256'].shape == (n, 32), 'Variant artifact shape mismatch.')
    hashes = np.ascontiguousarray(data['variant_sha256']).view('V32').ravel()
    need(len(np.unique(hashes)) == n, 'Duplicate immutable variant hash.')
    need(data['X'].shape == (len(data['member_lengths']), len(data['names'])) and
         int(data['gene_lengths'].sum()) == len(data['X']) and
         np.all(data['gene_lengths'] > 0) and np.all(data['member_lengths'] > 0) and
         int(data['member_lengths'].sum()) == len(data['member_flat']),
         'Packed gene/team membership shape mismatch.')
    need(np.all((data['member_flat'] >= 0) & (data['member_flat'] < n)) and
         len(np.unique(data['member_flat'])) == n, 'Invalid assigned membership indices.')
    data['sample_sha256'] = manifest['sample_sha256']
    data['manifest_sha256'] = done['manifest_sha256']
    data['boundary_diagnostics'] = diagnostics['boundary_diagnostics']
    return data

def descriptive_correlations(grouped, names, first):
    from scipy.stats import rankdata
    result = []
    for j, candidate in enumerate(first):
        if not candidate['candidate']:
            continue
        values = []
        for x, y in grouped:
            ok = np.isfinite(x[:, j]) & np.isfinite(y)
            if int(ok.sum()) < MIN_TEAMS:
                continue
            a, b = rankdata(x[ok, j]), rankdata(y[ok])
            a, b = a - a.mean(), b - b.mean()
            scale = np.linalg.norm(a) * np.linalg.norm(b)
            if scale > 0:
                values.append(float(a.dot(b) / scale))
        result.append(dict(annotation=names[j], eligible_genes=len(values),
                           observed=float(np.mean(values)) if values else None))
    return result

def run_pool(chrs):
    need(sorted(chrs, key=int) == ['2', '12'], 'Main gate requires exactly pool 2 12.')
    chrs = ['2', '12']
    out = BASE / 'out/pooled_2_12'
    lock = fresh_directory(out)
    provenance = implementation()
    built = [read_built(ch, provenance) for ch in chrs]
    need(built[0]['sample_sha256'] == built[1]['sample_sha256'],
         'Chromosomes used different LD sample sets or order.')
    names = sorted(set().union(*(set(map(str, b['names'])) for b in built)))
    all_grouped, chromosome_groups, flat_members = [], [], []
    variant_offset = 0
    for b in built:
        old = list(map(str, b['names']))
        need(len(old) == len(set(old)), 'Duplicate annotation names.')
        x = np.zeros((len(b['X']), len(names)), dtype=np.float64)
        for j, name in enumerate(names):
            if name in old:
                x[:, j] = b['X'][:, old.index(name)]
            else:
                need(name.startswith(('ccre_class=', 'rep_class=')),
                     'Missing noncategorical annotation column across chromosomes.')
        starts = np.r_[0, np.cumsum(b['gene_lengths'])]
        local = [(x[starts[g]:starts[g + 1]], np.full(b['gene_lengths'][g], np.nan))
                 for g in range(len(b['gene_lengths']))]
        all_grouped.extend(local)
        chromosome_groups.append(local)
        offsets = np.r_[0, np.cumsum(b['member_lengths'])]
        perteam = [b['member_flat'][offsets[t]:offsets[t + 1]] + variant_offset
                   for t in range(len(b['member_lengths']))]
        flat_members.extend([perteam[starts[g]:starts[g + 1]]
                             for g in range(len(b['gene_lengths']))])
        variant_offset += len(b['maf'])
    first = measure_one(all_grouped, names)
    beta = np.concatenate([b['absbeta'] for b in built])
    maf = np.concatenate([b['maf'] for b in built])
    distance = np.concatenate([b['lead'] for b in built])
    target, regression = residual_target(beta, maf, distance)
    attach_team_targets(all_grouped, flat_members, target)
    del beta, maf, distance, target, flat_members
    second = measure_two(all_grouped, names, first)
    descriptive = {}
    for ch, groups in zip(chrs, chromosome_groups):
        m1 = measure_one(groups, names)
        for row in m1:
            row.pop('candidate')
        descriptive[ch] = dict(measurement_1=m1,
                              measurement_2=descriptive_correlations(groups, names, first),
                              inferential=False)
    passed = [r['annotation'] for r in second if r['passed']]
    decision = ('PROCEED_TO_PHI_TRAINING' if passed else
                'NO_DISCRIMINABLE_EFFECT_HEIGHT_CHR2+12')
    report = dict(provenance=provenance, chromosomes=[2, 12],
                  main_scope='pooled supplied universe within protein_coding body +/-3kb',
                  input_manifest_sha256=[b['manifest_sha256'] for b in built],
                  sample_sha256=built[0]['sample_sha256'], regression=regression,
                  measurement_1=first, measurement_2=second,
                  chromosome_descriptive=descriptive, genes=len(all_grouped),
                  boundary_diagnostics={ch: b['boundary_diagnostics']
                                        for ch, b in zip(chrs, built)},
                  interpretation_limits=[
                      'Cross-block LD is ignored by the sealed team definition.',
                      'Boundary windows do not measure all omitted high-LD pairs.',
                      'Splitting dependent signals need not make label permutations conservative.',
                      'Pooled nominal permutation/BH decision is not a calibration proof.'],
                  assigned_unique_variants=variant_offset)
    final = dict(chromosomes=[2, 12], pooled=True, passed_annotations=passed,
                 decision=decision, candidates=len(second),
                 estimable_candidates=sum(r['status'] == 'tested' for r in second),
                 next_action='skip_genome_gate_and_train_stage_1' if passed else 'stop',
                 negative_interpretation=None if passed else
                 '주석으로 구별 가능한 효과 없음 (키, chr12+chr2)')
    need(implementation() == provenance and
         all(input_signature(ch) == json.loads((PRIVATE / ('chr%s' % ch) /
             'manifest.json').read_text())['inputs'] for ch in chrs),
         'Inputs changed during pooled calculation.')
    atomic_json(out / 'diagnostics.json', report)
    atomic_json(out / 'decision.json', final)
    lock.rmdir()
    atomic_json(out / 'gate.done', dict(status='complete', decision=decision,
                outputs={name: digest_file(out / name)
                         for name in ('diagnostics.json', 'decision.json')}))
    log('pool_complete', passed_annotations=len(passed))

def main():
    global np, LOCK_HANDLE
    need(len(sys.argv) >= 3 and sys.argv[1] in ('chr', 'pool'),
         'Usage: final_phi_gate_r3.py chr 2|12 OR pool 2 12.')
    mode, chrs = sys.argv[1], [chrom(x) for x in sys.argv[2:]]
    need((mode == 'chr' and len(chrs) == 1 and chrs[0] in ('2', '12')) or
         (mode == 'pool' and sorted(chrs, key=int) == ['2', '12']),
         'Only chr 2, chr 12, and pooled 2 12 are within the seal.')
    cap = limit_resources()
    with open(os.devnull, 'w') as null:
        os.dup2(null.fileno(), 2)
    guard()
    owned(BASE)
    lockpath = BASE / '.final_phi_gate.lock'
    need(not lockpath.is_symlink(), 'Global lock may not be a symlink.')
    fd = os.open(str(lockpath), os.O_CREAT | os.O_RDWR | os.O_NOFOLLOW, 0o600)
    need(os.fstat(fd).st_uid == os.getuid(), 'Global lock is not owned.')
    LOCK_HANDLE = os.fdopen(fd, 'r+')
    try:
        fcntl.flock(LOCK_HANDLE, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise GateError('Another final chr/pool job is already active.')
    make_owned_directory(PRIVATE)
    scratch = PRIVATE / 'scratch'
    make_owned_directory(scratch)
    for name in ('TMPDIR', 'TMP', 'TEMP'):
        os.environ[name] = str(scratch)
    os.nice(10)
    import numpy as numpy_module
    np = numpy_module
    monitor = threading.Thread(target=heartbeat, daemon=True)
    monitor.start()
    log('start', mode=mode, address_space_cap_gib=cap / GIB,
        process_count=1, affinity_cpus=len(os.sched_getaffinity(0)),
        blas_threads=1, dosage_readers=1)
    try:
        if mode == 'chr':
            run_chr(chrs[0])
        else:
            run_pool(chrs)
    finally:
        STOP.set()
        monitor.join(timeout=3)
        LOCK_HANDLE.close()

if __name__ == '__main__':
    try:
        main()
    except GateError as error:
        log('failed', reason=str(error))
        sys.exit(1)
    except (Exception, KeyboardInterrupt) as error:
        log('failed', reason='Unexpected error; incomplete outputs are not resumable.',
            exception_type=type(error).__name__)
        sys.exit(1)
