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
import subprocess
import sys
import time
import zipfile
from contextlib import contextmanager

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
FLANK_BP = 3000
E2G_SCORE = 0.6
VARIANCE_RATIO_MIN = 0.5
PERMUTATIONS = 1000
FDR_ALPHA = 0.10
THREADS = 4
GIB = 1024 ** 3
AS_CAP = 36 * GIB
WORKING_BUDGET = 24 * GIB
QUERY_BP = 1000000
DOT_BLOCK = 256
PARQUET_BATCH = 16384
PROGRESS_GENES = 50
ROOT = Path(_config_path('${PROJECT_ROOT}/work'))
REF = ROOT / 'ref'
BCF = Path('bcftools')
GTF = REF / 'deductive/gencode.nochr.gtf.gz'
GPN = REF / 'gpnmsa/scores.tsv.bgz'
BBJ = REF / 'bbj/hum0197.v3.BBJ.Hei.v1.zip'
BBJ_MEMBER = 'hum0197.v3.BBJ.Hei.v1/GWASsummary_Height_Japanese_SakaueKanai2020.auto.txt.gz'
LEADS = ROOT / 'prs/out/hei/lead_snps.tsv'
AUDIT = collections.Counter()
np = None

class GateError(Exception):
    pass

def need(condition, message):
    if not condition:
        raise GateError(message)

def log(stage, **counts):
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024
    record = dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  stage=stage, rss_peak_gib=round(rss / GIB, 3))
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
             'deeprvat_associate', 'deeprvat-associate')
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

@contextmanager
def command_lines(args, samples=None):
    proc = subprocess.Popen([str(x) for x in args],
                            stdin=subprocess.PIPE if samples is not None else subprocess.DEVNULL,
                            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                            text=True, bufsize=1 << 20)
    try:
        if samples is not None:
            proc.stdin.write('\n'.join(samples) + '\n')
            proc.stdin.close()
        yield proc.stdout
        need(proc.wait() == 0, 'Child process failed; raw stderr suppressed for privacy.')
    finally:
        proc.stdout.close()
        if proc.stdin is not None and not proc.stdin.closed:
            proc.stdin.close()
        if proc.poll() is None:
            proc.terminate()
            try:
                proc.wait(timeout=10)
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.wait()

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
    import pyarrow as pa
    import pyarrow.parquet as pq
    pa.set_cpu_count(1)
    pa.set_io_thread_count(1)
    log('assignment_start', variants=len(positions))
    genes = {}
    with gzip.open(GTF, 'rt') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            row = line.rstrip('\n').split('\t')
            if row[2] != 'gene':
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]*)"', row[8]))
            if attrs.get('gene_type') != 'protein_coding':
                continue
            c = chrom(row[0])
            if c not in set(map(str, range(1, 23))):
                continue
            gid = attrs['gene_id'].split('.')[0]
            need(gid not in genes, 'Duplicate version-stripped GTF gene ID.')
            a, b = int(row[3]), int(row[4])
            need(0 < a <= b, 'Invalid GTF interval.')
            genes[gid] = (c, a, b)
    intervals = {g: [(max(1, a - FLANK_BP), b + FLANK_BP)]
                 for g, (c, a, b) in genes.items() if c == ch}
    files = sorted((REF / 're2g_all/ot_e2g').glob('*.parquet'))
    need(files and intervals, 'No OT files or chromosome protein-coding genes.')
    columns = ['geneId', 'chromosome', 'start', 'end', 'score']
    seen_elements, mismatch_elements = set(), set()
    for file_no, path in enumerate(files, 1):
        with pq.ParquetFile(path) as pf:
            need(set(columns) <= set(pf.schema_arrow.names), 'OT schema mismatch.')
            for batch in pf.iter_batches(batch_size=PARQUET_BATCH, columns=columns, use_threads=False):
                data = batch.to_pydict()
                for g, c, a, b, score in zip(*(data[n] for n in columns)):
                    if score is None:
                        AUDIT['ot_null_score'] += 1
                        continue
                    score = float(score)
                    need(math.isfinite(score) and 0 <= score <= 1, 'Invalid OT score.')
                    if score < E2G_SCORE:
                        AUDIT['ot_below_score'] += 1
                        continue
                    need(isinstance(g, str), 'Invalid OT gene ID.')
                    g, c = g.split('.')[0], chrom(c)
                    need(type(a) is int and type(b) is int and 0 <= a < b, 'Invalid OT BED interval.')
                    if g not in genes:
                        AUDIT['ot_unmatched_gene_rows'] += 1
                        continue
                    need(c in set(map(str, range(1, 23))), 'Nonautosomal OT link to autosomal GTF gene.')
                    AUDIT['ot_gtf_joined_rows'] += 1
                    if c != genes[g][0]:
                        AUDIT['ot_gtf_chrom_mismatch_rows'] += 1
                        mismatch_elements.add((g, c, a, b))
                        continue
                    if c != ch:
                        continue
                    element = (g, a, b)
                    if element in seen_elements:
                        AUDIT['ot_duplicate_gene_element_rows'] += 1
                        continue
                    seen_elements.add(element)
                    intervals[g].append((a + 1, b))
        log('assignment_ot_file', files_done=file_no, files_total=len(files), elements=len(seen_elements))
    AUDIT['ot_gtf_chrom_mismatch_elements'] = len(mismatch_elements)
    need(AUDIT['ot_gtf_chrom_mismatch_rows'] * 1000 <= AUDIT['ot_gtf_joined_rows'],
         'OT/GTF chromosome mismatch exceeds the inherited 0.1 percent gate.')
    need(seen_elements, 'Zero usable chromosome OT/GTF join.')
    order = np.argsort(positions, kind='stable')
    sorted_pos = positions[order]
    assigned = []
    for gid in sorted(intervals):
        merged = []
        for a, b in sorted(intervals[gid]):
            if merged and a <= merged[-1][1] + 1:
                merged[-1] = (merged[-1][0], max(b, merged[-1][1]))
            else:
                merged.append((a, b))
        pieces = [order[np.searchsorted(sorted_pos, a, side='left'):
                        np.searchsorted(sorted_pos, b, side='right')] for a, b in merged]
        members = np.sort(np.concatenate(pieces))
        if len(members):
            need(len(members) == len(np.unique(members)), 'Duplicate within-gene variant assignment.')
            assigned.append(members)
    need(assigned, 'No gene assignments.')
    AUDIT['genes_assigned'] = len(assigned)
    AUDIT['gene_variant_instances'] = sum(len(x) for x in assigned)
    log('assignment_done', genes=len(assigned), instances=AUDIT['gene_variant_instances'])
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
                candidates.append(p0 - intervals[i][1] + 1)
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
    for tf in top:
        columns['tf_' + tf] = np.zeros(len(keys), dtype=np.float64)
    columns['tf_n'] = np.zeros(len(keys), dtype=np.float64)
    for n, p in enumerate(positions):
        p0 = int(p) - 1
        j = bisect.bisect_right(starts, p0) - 1
        seen = set()
        while j >= 0 and p0 - starts[j] < 20000:
            a, b, tf = iv[j]
            if a <= p0 < b:
                seen.add(tf)
            j -= 1
        columns['tf_n'][n] = len(seen)
        for tf in seen.intersection(top):
            columns['tf_' + tf][n] = 1
    del iv, starts
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
            records = {}
            for line in tab.fetch(contigs[0], min(bypos) - 1, max(bypos)):
                row = line.split('\t')
                need(len(row) == 5 and chrom(row[0]) == ch, 'GPN five-column schema mismatch.')
                p = int(row[1])
                if p not in bypos:
                    continue
                k = (p, row[2].upper(), row[3].upper())
                need(k not in records, 'Duplicate GPN allele key.')
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
            if block_no % 10 == 0:
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

def residual_target(beta, mafs, positions, leads):
    log('residuals_start', variants=len(beta))
    lp = np.array(sorted(set(leads.values())), dtype=np.int64)
    right = np.searchsorted(lp, positions)
    distance = np.minimum(np.abs(positions - lp[np.maximum(0, right - 1)]),
                          np.abs(positions - lp[np.minimum(len(lp) - 1, right)]))
    valid = np.isfinite(beta) & (beta > 0) & (mafs > 0)
    AUDIT['bbj_zero_beta'] = int(np.sum(beta == 0))
    AUDIT['matched_nonpositive_maf'] = int(np.sum(np.isfinite(beta) & (mafs <= 0)))
    need(int(valid.sum()) > 1 + 2 * (INTERNAL_KNOTS + 1), 'Insufficient valid BBJ targets for spline OLS.')
    x = np.column_stack([np.ones(int(valid.sum())), natural_spline(np.log10(mafs[valid])),
                         natural_spline(np.log10(distance[valid] + 1.0))])
    y = np.log(beta[valid])
    coef, _, rank, singular = np.linalg.lstsq(x, y, rcond=None)
    need(rank == x.shape[1], 'Spline design is rank deficient.')
    residual = y - x.dot(coef)
    need(np.isfinite(residual).all() and np.var(residual) > 0, 'Invalid or constant residual target.')
    result = np.full(len(beta), np.nan)
    result[valid] = residual
    AUDIT['residual_eligible_variants'] = int(valid.sum())
    log('residuals_done', fitted=int(valid.sum()), design_columns=int(x.shape[1]), design_rank=int(rank))
    return result, dict(n_fitted=int(valid.sum()), design_rank=int(rank),
                        design_columns=int(x.shape[1]), condition_number=float(singular[0] / singular[-1]),
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
    indices = sorted(random.Random(SEED).sample(range(len(samples)), SAMPLE_N))
    chosen = [samples[i] for i in indices]
    log('dosage_header_done', available_samples=len(samples), selected_samples=len(chosen),
        normalized_contig_matches=len(contigs))
    return path, contigs[0], chosen, ds_number

def gene_dosage(ids, keys, path, contig, samples, ds_number):
    n = len(ids)
    projected = 3 * n * SAMPLE_N * 8 + 4 * (n * (n - 1) // 2) * 8
    need(projected <= WORKING_BUDGET, 'Gene exact-LD working set exceeds memory budget; no truncation allowed.')
    matrix = np.empty((n, SAMPLE_N), dtype=np.float64)
    local = {keys[int(i)]: j for j, i in enumerate(ids)}
    found = np.zeros(n, dtype=bool)
    regions = collections.defaultdict(list)
    for k in local:
        regions[(k[1] - 1) // QUERY_BP].append(k[1])
    for positions in regions.values():
        region = '%s:%d-%d' % (contig, min(positions), max(positions))
        args = [BCF, 'query', '-S', '-', '-r', region, '--regions-overlap', '0',
                '-f', '%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n', path]
        with command_lines(args, samples=samples) as lines:
            for line in lines:
                row = line.rstrip('\n').split('\t')
                need(len(row) == 4 + SAMPLE_N, 'Dosage sample cardinality changed.')
                c, p, r = chrom(row[0]), int(row[1]), row[2].upper()
                alts = row[3].upper().split(',')
                wanted = [(ai, local[(c, p, r, a)]) for ai, a in enumerate(alts)
                          if (c, p, r, a) in local]
                if not wanted:
                    continue
                need(len(alts) == 1 or ds_number == 'A', 'Multiallelic DS lacks Number=A semantics.')
                split_values = [v.split(',') for v in row[4:]]
                need(all(len(v) == len(alts) or v == ['.'] for v in split_values),
                     'DS vector length differs from allele count.')
                for ai, j in wanted:
                    need(not found[j], 'Duplicate dosage allele key in queried regions.')
                    values = np.array([float(v[ai]) if v != ['.'] and v[ai] != '.' else np.nan
                                       for v in split_values], dtype=np.float64)
                    good = np.isfinite(values)
                    need(good.any() and np.all((values[good] >= 0) & (values[good] <= 2)),
                         'All-missing or invalid dosage.')
                    AUDIT['ds_missing_cells_gene_instances'] += int((~good).sum())
                    values[~good] = values[good].mean()
                    matrix[j] = values
                    found[j] = True
    need(found.all(), 'Incomplete dosage join; cannot form exact teams.')
    matrix -= matrix.mean(axis=1, keepdims=True)
    norms = np.sqrt(np.einsum('ij,ij->i', matrix, matrix))
    variable = norms > 0
    AUDIT['ds_constant_gene_instances'] += int((~variable).sum())
    matrix[variable] /= norms[variable, None]
    matrix[~variable] = 0
    return matrix, variable

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
            d = np.where(r2 > LD_R2, 1.0 - r2, np.maximum(1.0 - r2, 1.0 - LD_R2))
            offset = n * j - j * (j + 1) // 2
            dist[offset:offset + len(d)] = d
    tree = linkage(dist, method='complete')
    labels = fcluster(tree, t=np.nextafter(1.0 - LD_R2, -np.inf), criterion='distance')
    groups = collections.defaultdict(list)
    for i, label in enumerate(labels):
        groups[int(label)].append(i)
    teams = [np.array(g, dtype=np.int64) for g in sorted(groups.values(), key=lambda g: g[0])]
    for members in teams:
        for i, j in enumerate(members[:-1]):
            offsets = n * j - j * (j + 1) // 2 + members[i + 1:] - j - 1
            need(np.all(dist[offsets] < 1.0 - LD_R2), 'Complete-link clique invariant failed.')
    return teams

def build_teams(assignments, keys, features, dosage):
    log('teams_start', genes=len(assignments))
    grouped = []
    memberships = []
    sizes = collections.Counter()
    for g, ids in enumerate(assignments, 1):
        matrix, variable = gene_dosage(ids, keys, *dosage)
        teams = complete_teams(matrix)
        del matrix, variable
        x = np.full((len(teams), features.shape[1]), np.nan)
        y = np.full(len(teams), np.nan)
        gene_memberships = []
        for t, members in enumerate(teams):
            variants = ids[members]
            gene_memberships.append(variants)
            sub = features[variants]
            valid = np.isfinite(sub)
            count = valid.sum(axis=0)
            x[t] = np.divide(np.where(valid, sub, 0).sum(axis=0), count,
                             out=np.full(features.shape[1], np.nan), where=count > 0)
            sizes[len(members)] += 1
        grouped.append((x, y))
        memberships.append(gene_memberships)
        if g % PROGRESS_GENES == 0 or g == len(assignments):
            log('teams_progress', genes_done=g, genes_total=len(assignments), teams=sum(sizes.values()))
            guard()
    need(sum(sizes.values()) > 0, 'No teams.')
    AUDIT['teams'] = sum(sizes.values())
    log('teams_done', genes=len(grouped), teams=AUDIT['teams'])
    return grouped, dict(team_count=sum(sizes.values()), singleton_count=sizes[1],
                         largest_team=max(sizes), mean_team_size=sum(k * v for k, v in sizes.items()) / sum(sizes.values())), memberships

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
            for b in range(PERMUTATIONS):
                null_sum[b] += float(xr.dot(rng.permutation(yr)))
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

def main():
    global np
    need(len(sys.argv) == 2 and chrom(sys.argv[1]) in set(map(str, range(1, 23))),
         'Usage: draft_b_phi_gate.py <autosomal chromosome 1..22>.')
    ch = chrom(sys.argv[1])
    out = ROOT / ('phi_gate/out/chr%s' % ch)
    need(not out.is_symlink(), 'Output root may not be a symlink.')
    if out.exists():
        need(not any(p.name.endswith('.done') or p.name == 'done' for p in out.iterdir()),
             'Done marker exists; refusing to rerun.')
        need(not any(out.iterdir()), 'Partial/nonempty output directory; manual review required.')
    cap = limit_resources()
    guard()
    owned(ROOT / 'phi_gate')
    out.parent.mkdir(exist_ok=True)
    owned(out.parent)
    out.mkdir(exist_ok=True)
    owned(out)
    lock = out / '.running'
    lock.mkdir()
    for name in ('TMPDIR', 'TMP', 'TEMP'):
        os.environ[name] = str(out)
    with open(os.devnull, 'w') as null:
        os.dup2(null.fileno(), 2)
    import numpy as numpy_module
    np = numpy_module
    log('start', chromosome=int(ch), address_space_cap_gib=cap / GIB)
    keys, positions, mafs, index = read_universe(ch)
    leads = read_leads(ch)
    lifted = read_lifts(ch, keys, positions, index, leads)
    assignments = assign_genes(ch, positions)
    names, features, coverage = annotations(ch, keys, positions, lifted)
    dosage = dosage_header(ch)
    grouped, team_summary, memberships = build_teams(assignments, keys, features, dosage)
    del features, lifted, dosage
    first = measure_one(grouped, names)
    beta, matched = match_bbj(ch, keys, index)
    target, regression = residual_target(beta, mafs, positions, leads)
    attach_team_targets(grouped, memberships, target)
    del beta, memberships
    second = measure_two(grouped, names, first)
    passed = [r['annotation'] for r in second if r['passed']]
    decision = 'PROCEED_TO_PHI_TRAINING' if passed else 'NO_DISCRIMINABLE_EFFECT_HEIGHT_CHR%s' % ch
    report = dict(schema_version=1, implementation='draft_b_r1', chromosome=int(ch),
                  scope='supplied filtered universe on this chromosome only',
                  audit=dict(AUDIT), bbj_matching_rate=float(matched.mean()),
                  residual_eligible_rate=AUDIT['residual_eligible_variants'] / len(keys),
                  annotations=coverage, regression=regression, teams=team_summary,
                  measurement_1=first, measurement_2=second,
                  defaults=dict(samples=SAMPLE_N, seed=SEED, min_teams=MIN_TEAMS,
                                internal_knots=INTERNAL_KNOTS, knot_quantiles=KNOT_QUANTILES,
                                palindromic_policy=PALINDROMIC_POLICY, ds_missing=DS_MISSING_POLICY,
                                feature_missing=FEATURE_MISSING_POLICY),
                  sealed=dict(qc_r2=R2_QC, ld_r2_strict=LD_R2, permutations=PERMUTATIONS,
                              variance_ratio_min=VARIANCE_RATIO_MIN, fdr_strict=FDR_ALPHA))
    final = dict(chromosome=int(ch), passed_annotations=passed, decision=decision,
                 candidates=len(second), estimable_candidates=sum(r['status'] == 'tested' for r in second),
                 bbj_matching_rate=float(matched.mean()), report='diagnostics.json')
    atomic_json(out / 'diagnostics.json', report)
    atomic_json(out / 'decision.json', final)
    lock.rmdir()
    atomic_json(out / 'phi_gate.done', dict(status='complete', decision=decision,
                                           outputs=['diagnostics.json', 'decision.json']))
    log('complete', chromosome=int(ch), passed_annotations=len(passed))

if __name__ == '__main__':
    try:
        main()
    except GateError as error:
        log('failed', reason=str(error))
        sys.exit(1)
    except (Exception, KeyboardInterrupt) as error:
        log('failed', reason='Unexpected error; no done marker should be consumed without outputs.',
            exception_type=type(error).__name__)
        sys.exit(1)
