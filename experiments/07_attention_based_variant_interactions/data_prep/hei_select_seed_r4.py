#!/usr/bin/env python
import os as _os
if not _os.environ.get("PROJECT_ROOT", "").strip():
    raise ValueError("PROJECT_ROOT must be set and nonblank")
import bisect, collections, gzip, json, math, os, re, resource, shutil, subprocess, sys, time
from contextlib import contextmanager
from pathlib import Path

W = Path((_os.environ["PROJECT_ROOT"] + '/work/prs')); R = Path((_os.environ["PROJECT_ROOT"] + '/work/ref'))
C = W / 'out/hei/clump'; S = W / 'out/hei/sel_seed'; TMP = Path((_os.environ["PROJECT_ROOT"] + '/work/tmp'))
B = Path(_os.environ.get("BCFTOOLS_BIN", "bcftools"))
GTF = R / 'deductive/gencode.nochr.gtf.gz'; OT = R / 're2g_all/ot_e2g'
CHS = tuple(map(str, range(1, 23))); CAP = 256; GIB = 1024 ** 3
SEEDS = ('chondrocyte', 'osteoblast', 'osteocyte', 'femur', 'mesenchymal stem cell',
         'dedifferentiated amniotic fluid mesenchymal stem cell', 'stromal cell of bone marrow',
         'HS-5', 'HS-27A', 'limb', 'left forelimb', 'right forelimb', 'left hindlimb', 'right hindlimb')
class SelectionError(Exception): pass
def need(ok, message):
    if not ok: raise SelectionError(message)
def chrom(value):
    value = str(value).strip(); value = value[3:] if value.lower().startswith('chr') else value
    return str(int(value)) if value.isdigit() else value
def key19(value):
    a = value.strip().split(':')
    need(len(a) == 4 and a[1].isdigit() and int(a[1]) > 0, 'Invalid hg19 key.')
    need(all(x and x != '.' and ';' not in x for x in a[2:]), 'Invalid allele key.')
    return ':'.join((chrom(a[0]), str(int(a[1])), a[2], a[3]))
def number(value, low, high):
    value = float(value); need(math.isfinite(value) and low <= value <= high, 'Numeric contract failed.')
    return value
def owned(path):
    need(path.is_dir() and path.stat().st_uid == os.getuid() and
         str(path.resolve()).startswith((_os.environ["PROJECT_ROOT"] + '/')), 'Owned /data directory required.')
def limits():
    owned(W); owned(S.parent); owned(TMP); sys.dont_write_bytecode = True; os.umask(0o077)
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    cap = min([120 * 1000**3] + [v for v in (soft, hard) if v != resource.RLIM_INFINITY])
    resource.setrlimit(resource.RLIMIT_AS, (cap, cap))
    os.sched_setaffinity(0, set(sorted(os.sched_getaffinity(0))[:4]))
    for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS',
              'VECLIB_MAXIMUM_THREADS', 'ARROW_NUM_THREADS'): os.environ[k] = '1'
    for k in ('TMPDIR', 'TMP', 'TEMP'): os.environ[k] = str(TMP)
    return cap / GIB
def snapshot():
    records = {}; heavy = ('deepsea', 'deepripe', 'sparse_genotype', 'sparse-genotype',
                          'step2_spa', 'saige_step2', 'saige-step2', 'step2_tests', 'deeprvat_associate', 'deeprvat-associate')
    for p in Path('/proc').iterdir():
        if not p.name.isdigit() or int(p.name) == os.getpid(): continue
        try:
            uid = p.stat().st_uid; raw = (p / 'stat').read_text(); a = raw[raw.rfind(')') + 2:].split()
            records[(p.name, a[19])] = (uid, int(a[11]) + int(a[12]), max(0, int(a[21])) * os.sysconf('SC_PAGE_SIZE'))
            tokens = [os.path.basename(t).lower() for t in (p / 'cmdline').read_text().split('\0') if t]
            busy = any(n in t for n in heavy for t in tokens) or (
                any('deeprvat' in t for t in tokens) and any(t in ('associate', 'association') for t in tokens))
            need(not busy, 'Resource guard: heavyweight stage active.')
        except (FileNotFoundError, ProcessLookupError): continue
    return records
def guard():
    before = snapshot(); started = time.monotonic(); time.sleep(1); after = snapshot()
    elapsed = time.monotonic() - started; cpu = collections.Counter(); ram = collections.Counter()
    for identity, (uid, ticks, rss) in after.items():
        if uid == os.getuid(): continue
        ram[uid] += rss
        cpu[uid] += max(0, ticks - before.get(identity, (uid, 0, 0))[1]) / os.sysconf('SC_CLK_TCK') / elapsed * 100
    need(max(ram.values(), default=0) <= 32 * GIB and sum(ram.values()) <= 64 * GIB, 'Other-user RAM limit.')
    need(max(cpu.values(), default=0) <= 400 and sum(cpu.values()) <= 800, 'Other-user CPU limit.')
    mem = dict(line.split(':', 1) for line in Path('/proc/meminfo').read_text().splitlines())
    need(int(mem['MemAvailable'].split()[0]) * 1024 >= 80 * GIB, 'Available RAM below 80 GiB.')
    need(shutil.disk_usage(W).free >= 10 * GIB, 'Output disk below 10 GiB free.')
    if Path('/dev/nvidia0').exists():
        proc = subprocess.run(['nvidia-smi', '--query-compute-apps=pid', '--format=csv,noheader,nounits'],
                              stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, timeout=15)
        need(proc.returncode == 0, 'GPU inspection failed.')
        for pid in proc.stdout.split():
            need(pid.isdigit(), 'GPU process schema failed.')
            try: need((Path('/proc') / pid).stat().st_uid == os.getuid(), 'Other-user GPU is busy.')
            except FileNotFoundError: pass
def rows(path, required):
    with open(path) as f:
        h = f.readline().lstrip('#').split()
        need(len(h) == len(set(h)) and set(required.split()) <= set(h), 'Cache header failed.')
        for line in f:
            if not line.strip(): continue
            a = line.split(); need(len(a) == len(h), 'Cache row failed.'); yield dict(zip(h, a))
def caches():
    leads = {}; proxies = {}; audit = collections.Counter(cs_rows_non_autosomal_skipped=0); groups = collections.defaultdict(dict)
    with open(C / 'split.done') as f: nsig = json.load(f)['n_sig_nonMHC']
    need(type(nsig) is int and nsig >= 0, 'Split completion cache failed.')
    for ch in CHS:
        local = {}
        for r in rows(C / f'chr{ch}.clumps', 'CHROM POS ID P'):
            k = key19(r['ID']); p = int(k.split(':')[1])
            need(chrom(r['CHROM']) == ch == k.split(':')[0] and int(r['POS']) == p, 'Lead coordinate mismatch.')
            need(k not in local and not (ch == '6' and 25000000 <= p <= 34000000), 'Duplicate or MHC lead.')
            local[k] = number(r['P'], 0, 5e-8)
        leads[ch] = local; px = collections.defaultdict(set)
        if local:
            for r in rows(C / f'ld_chr{ch}.vcor', 'ID_A ID_B UNPHASED_R2'):
                if number(r['UNPHASED_R2'], 0, 1) < 0.5: audit['ld_rows_below_threshold'] += 1; continue
                a, b = key19(r['ID_A']), key19(r['ID_B'])
                need(a.split(':')[0] == ch == b.split(':')[0], 'LD chromosome mismatch.')
                for lead, proxy in ((a, b), (b, a)):
                    if lead in local and proxy != lead:
                        audit['duplicate_ld_edges'] += int(proxy in px[lead]); px[lead].add(proxy)
        proxies[ch] = px
    need(any(leads.values()), 'No independent leads.')
    for r in rows(C / 'susie_cs.tsv', 'chr pos19 a1 a2 region cs_id pip'):
        ch = chrom(r['chr'])
        if ch == 'X': audit['cs_rows_non_autosomal_skipped'] += 1; continue
        need(ch in CHS, 'Invalid CS chromosome.')
        k = key19(f"{ch}:{r['pos19']}:{r['a1']}:{r['a2']}")
        need(r['cs_id'] != '-1', 'Invalid CS group.')
        group = groups[(ch, r['region'], r['cs_id'])]; need(k not in group, 'Duplicate CS variant.')
        group[k] = number(r['pip'], 0, 1)
    assigned = {ch: collections.defaultdict(set) for ch in CHS}; nassigned = 0
    for (ch, region, cid), vs in groups.items():
        if not leads[ch]: continue
        top = int(max(vs, key=vs.get).split(':')[1])
        d, lead = min((abs(top - int(k.split(':')[1])), k) for k in leads[ch])
        if d <= 500000:
            nassigned += 1
            for k in vs:
                c, p, a, b = k.split(':'); assigned[ch][lead].update((k, f'{c}:{p}:{b}:{a}'))
    return leads, proxies, assigned, dict(n_cs=len(groups), n_cs_assigned=nassigned), nsig, audit
class Intervals:
    def __init__(self, pairs):
        merged = []
        for a, b in sorted(pairs):
            if merged and a <= merged[-1][1] + 1: merged[-1] = (merged[-1][0], max(b, merged[-1][1]))
            else: merged.append((a, b))
        self.starts = [a for a, b in merged]; self.ends = [b for a, b in merged]
    def contains(self, p):
        i = bisect.bisect_right(self.starts, p) - 1; return i >= 0 and p <= self.ends[i]
def annotation():
    genes = {}; cds = {ch: [] for ch in CHS}
    with gzip.open(GTF, 'rt') as f:
        for line in f:
            if line.startswith('#'): continue
            a = line.rstrip('\n').split('\t'); need(len(a) == 9, 'Invalid GTF row.'); ch = chrom(a[0])
            if ch not in CHS or a[2] not in ('gene', 'CDS'): continue
            s, e = int(a[3]), int(a[4]); need(1 <= s <= e, 'Invalid GTF interval.')
            if a[2] == 'CDS': cds[ch].append((s, e)); continue
            attrs = dict(re.findall(r'(\w+) "([^"]*)"', a[8]))
            if attrs.get('gene_type') != 'protein_coding': continue
            gid = attrs['gene_id'].split('.')[0]; name = attrs['gene_name']
            need(gid and name and gid not in genes, 'GTF gene ID is empty or duplicated.')
            genes[gid] = (ch, name, s, e)
    need(genes, 'No protein-coding genes.'); return genes, {ch: Intervals(v) for ch, v in cds.items()}
def elements(genes, files, tissue, audit):
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    pa.set_cpu_count(1); pa.set_io_thread_count(1)
    cols = ['biosampleName', 'geneId', 'chromosome', 'start', 'end', 'score']; best = {}; names = pa.array(SEEDS)
    mismatch_names = set()
    audit.update(ot_gtf_joined_rows=0, ot_gtf_chrom_mismatch_rows=0, ot_gtf_chrom_mismatch_elements=0)
    for path in files:
        with pq.ParquetFile(path) as pf:
            need(set(cols) <= set(pf.schema_arrow.names), 'Missing OT columns.')
            for batch in pf.iter_batches(batch_size=32768, columns=cols, use_threads=False):
                batch = batch.filter(pc.is_in(batch.column(batch.schema.get_field_index('biosampleName')), value_set=names))
                data = batch.to_pydict()
                for t, g, ch, a, b, score in zip(*(data[c] for c in cols)):
                    tissue[t]['raw_rows'] += 1
                    if score is None: tissue[t]['null_score_rows'] += 1; audit['null_score_rows'] += 1; continue
                    score = number(score, 0.6, 1)
                    need(isinstance(g, str) and g.split('.')[0], 'Invalid OT gene ID.'); g = g.split('.')[0]; ch = chrom(ch)
                    need(type(a) is int and type(b) is int and 0 <= a < b, 'Invalid OT BED interval.')
                    if g in genes:
                        need(ch in CHS, 'Invalid or non-autosomal OT chromosome for autosomal GTF gene.')
                        audit['ot_gtf_joined_rows'] += 1
                        if genes[g][0] != ch:
                            audit['ot_gtf_chrom_mismatch_rows'] += 1
                            mismatch_names.add(genes[g][1])
                    k = (ch, g, a, b); old = best.get(k, (0, set()))
                    audit['duplicate_gene_element_rows'] += int(k in best)
                    best[k] = (max(old[0], score), old[1] | {t})
    audit['ot_gtf_chrom_mismatch_genes'] = sorted(mismatch_names)
    need(audit['ot_gtf_chrom_mismatch_rows'] * 1000 <= audit['ot_gtf_joined_rows'],
         f"OT/GTF chromosome mismatch exceeds 0.1% of eligible joined rows ({audit['ot_gtf_chrom_mismatch_rows']}/{audit['ot_gtf_joined_rows']}).")
    el = collections.defaultdict(list)
    for (ch, g, a, b), (score, ts) in sorted(best.items()):
        for t in ts: tissue[t]['unique_gene_elements'] += 1
        if g not in genes: audit['unmatched_or_nonprotein_coding_elements'] += 1; continue
        if genes[g][0] != ch:
            audit['ot_gtf_chrom_mismatch_elements'] += 1; continue
        el[g].append((a, b, score, ts))
        for t in ts: tissue[t]['usable_gene_elements'] += 1
    audit['unique_seed_gene_elements'] = len(best); audit['usable_seed_gene_elements'] = sum(map(len, el.values()))
    need(el, 'OT/GTF join has zero usable elements.'); return el
def variants(ch, qc=False):
    path = R / f'lift38_keyed/chr{ch}.keyed38.vcf.gz'; proc = None
    if qc:
        proc = subprocess.Popen([str(B), 'query', '-i', 'INFO/R2>=0.3 && INFO/MAF>=0.0008',
                                 '-f', '%CHROM\t%POS\t%ID\n', str(path)], stdin=subprocess.DEVNULL,
                                stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
        f = proc.stdout
    else: f = gzip.open(path, 'rt')
    try:
        for line in f:
            if line.startswith('#'): continue
            a = line.rstrip('\n').split('\t', 3); need(len(a) == (3 if qc else 4), 'Invalid VCF query row.')
            if chrom(a[0]) != ch: continue
            k = key19(a[2]); p = int(a[1]); need(p > 0, 'Invalid hg38 position.')
            if k.split(':')[0] == ch: yield k, p
        if proc: need(proc.wait() == 0, 'bcftools failed.')
    finally:
        f.close()
        if proc:
            if proc.poll() is None:
                proc.terminate()
                try: proc.wait(timeout=5)
                except subprocess.TimeoutExpired: proc.kill()
            proc.wait()
def put(mapping, k, p, audit, label):
    need(k not in mapping or mapping[k] == p, 'Conflicting hg38 mapping.'); audit[label] += int(k in mapping); mapping[k] = p
def distance(p, s, e): return max(s - p, 0, p - e)
def process(ch, leads, px, cs, genes, el, cds, tissue, used_names):
    audit = collections.Counter({k: 0 for k in ('indep_leads', 'leads_not_lifted', 'gene_links', 'selected_genes',
                                'dist_shortage_leads', 'dist_shortage_tokens', 'dist_expected_tokens', 'dist_actual_tokens')})
    linked = {}; reasons = collections.Counter(); combos = collections.Counter()
    gr = []; dr = []; lr = []; unc = []; sources = []; hit_elements = {t: set() for t in SEEDS}
    if not leads: return gr, dr, lr, unc, sources, audit, reasons, combos
    gids = [g for g in genes if genes[g][0] == ch]; need(gids, 'No genes on lead chromosome.')
    wanted = set(leads).union(*(px.get(k, set()) | cs.get(k, set()) for k in leads)); lifted = {}
    stream = variants(ch)
    try:
        for k, p in stream:
            if k in wanted: put(lifted, k, p, audit, 'duplicate_link_lift_rows')
    finally: stream.close()
    audit.update(link_keys_requested=len(wanted), link_keys_lifted=len(lifted))
    for k, pv in leads.items():
        if k not in lifted: audit['leads_not_lifted'] += 1; continue
        p = lifted[k]; lv = {lifted[q] for q in ({k} | px.get(k, set()) | cs.get(k, set())) if q in lifted}
        sources.append((bool(px.get(k)), bool(cs.get(k)), len(lv))); lr.append((ch, k.split(':')[1], pv, k, p))
        nearest = min(gids, key=lambda g: distance(p, *genes[g][2:])); linked[k] = set()
        for g in gids:
            _, name, s, e = genes[g]; body = any(s - 3000 <= q <= e + 3000 for q in lv); ts = set()
            for a, b, score, tags in el.get(g, ()):
                if any(a < q <= b for q in lv):
                    ts.update(tags)
                    for t in tags: hit_elements[t].add((g, a, b))
            labels = [label for label, hit in (('body', body), ('element', bool(ts)), ('nearest', g == nearest)) if hit]
            if not labels: continue
            linked[k].add(g); reasons.update(labels); combos['+'.join(labels)] += 1
            for t in ts:
                tissue[t]['link_pairs'] += 1
                tissue[t]['element_only_link_pairs'] += int(not body and g != nearest)
                tissue[t]['exclusive_element_only_link_pairs'] += int(not body and g != nearest and len(ts) == 1)
    need(linked, 'No lead lifts on a populated chromosome.')
    for t in SEEDS: tissue[t]['linked_gene_elements'] += len(hit_elements[t])
    selected = sorted(set().union(*linked.values()), key=lambda g: (genes[g][1], g))
    names = [genes[g][1] for g in selected]
    need(len(names) == len(set(names)) and not (set(names) & used_names), 'Selected gene_name collision.'); used_names.update(names)
    iv = {g: [(max(1, genes[g][2] - 3000), genes[g][3] + 3000, 0.0)] +
          [(a + 1, b, sc) for a, b, sc, ts in el.get(g, ())] for g in selected}
    keep = Intervals([(a, b) for v in iv.values() for a, b, sc in v] +
                     [(max(1, lifted[k] - 1000000), lifted[k] + 1000000) for k in linked])
    candidates = {}; stream = variants(ch, True)
    try:
        for k, p in stream:
            if keep.contains(p) and not cds.contains(p): put(candidates, k, p, audit, 'duplicate_candidate_rows')
    finally: stream.close()
    items = sorted(candidates, key=lambda k: (candidates[k], k)); pos = [candidates[k] for k in items]; tokens = {}
    audit['unique_candidate_keys'] = len(candidates); del candidates
    for g in selected:
        _, name, s, e = genes[g]; scores = {}
        for a, b, sc in iv[g]:
            for i in range(bisect.bisect_left(pos, a), bisect.bisect_right(pos, b)): scores[i] = max(scores.get(i, 0), sc)
        unc.append(len(scores)); order = sorted(scores, key=lambda i: (-scores[i], distance(pos[i], s, e), pos[i], items[i]))[:CAP]
        tokens[g] = {items[i] for i in order}
        for rank, i in enumerate(order):
            direct = [sc for a, b, sc, ts in el.get(g, ()) if a < pos[i] <= b]
            need((direct or max(1, s - 3000) <= pos[i] <= e + 3000) and scores[i] == max(direct or [0]), 'Token membership/score mismatch.')
            gr.append((name, ch, items[i], pos[i], scores[i], rank))
    for k in linked:
        n = len(set().union(*(tokens[g] for g in linked[k]))); p = lifted[k]
        order = sorted(range(bisect.bisect_left(pos, max(1, p - 1000000)), bisect.bisect_right(pos, p + 1000000)),
                       key=lambda i: (abs(pos[i] - p), pos[i], items[i]))
        audit['dist_shortage_leads'] += int(len(order) < n); audit['dist_shortage_tokens'] += max(0, n - len(order))
        audit['dist_expected_tokens'] += n; audit['zero_token_leads'] += int(n == 0)
        dr.extend((k, items[i], pos[i]) for i in order[:n])
    audit.update(indep_leads=len(leads), gene_links=sum(map(len, linked.values())), selected_genes=len(selected))
    audit['dist_actual_tokens'] = len(dr)
    need(len(dr) + audit['dist_shortage_tokens'] == audit['dist_expected_tokens'], 'DIST count invariant failed.')
    need(sum(combos.values()) == audit['gene_links'], 'Link reason count invariant failed.')
    return gr, dr, lr, unc, sources, audit, reasons, combos
@contextmanager
def atomic(name, header=None):
    tmp = S / (name + '.tmp')
    with (gzip.open(tmp, 'xt') if name.endswith('.gz') else open(tmp, 'x')) as f:
        if header is not None: f.write(header + '\n')
        yield f
    with open(tmp, 'rb') as f: os.fsync(f.fileno())
    os.replace(tmp, S / name)
def write_rows(f, values):
    for row in values: f.write('\t'.join(map(str, row)) + '\n')
def sync_directory():
    fd = os.open(S, os.O_RDONLY | os.O_DIRECTORY)
    try: os.fsync(fd)
    finally: os.close(fd)
def quantiles(values, qs):
    values = sorted(values) or [0]; out = []
    for q in qs:
        i = (len(values) - 1) * q / 100; a = int(i); b = math.ceil(i)
        out.append(float(values[a] + (values[b] - values[a]) * (i - a)))
    return out
def manifest(files): return {str(p): (p.stat().st_size, p.stat().st_mtime_ns) for p in files}
def main():
    need(not S.exists() and not S.is_symlink(), 'sel_seed exists; no overwrite.')
    cap = limits(); guard(); files = sorted(OT.glob('*.parquet')); need(files, 'No OT parquet input.')
    inputs = files + [GTF, B, Path(__file__), C / 'split.done', C / 'susie_cs.tsv'] + list(C.glob('*.clumps')) + list(C.glob('*.vcor'))
    inputs += [R / f'lift38_keyed/chr{ch}.keyed38.vcf.gz' for ch in CHS]; initial = manifest(inputs)
    leads, px, cs, cssum, nsig, cache_audit = caches(); genes, cds = annotation()
    fields = 'raw_rows null_score_rows unique_gene_elements usable_gene_elements linked_gene_elements link_pairs element_only_link_pairs exclusive_element_only_link_pairs'
    tissue = {t: collections.Counter({k: 0 for k in fields.split()}) for t in SEEDS}
    ea = collections.Counter(null_score_rows=0, duplicate_gene_element_rows=0)
    el = elements(genes, files, tissue, ea); guard(); S.mkdir(mode=0o700)
    with atomic('seed_elements.tsv.gz', 'gene_id\tchr\tstart\tend\tscore\ttissues') as f:
        write_rows(f, ((g, genes[g][0], a, b, sc, '|'.join(sorted(ts))) for g in sorted(el) for a, b, sc, ts in el[g]))
    totals = collections.Counter(); gc = collections.Counter()
    reasons = collections.Counter(body=0, element=0, nearest=0); combos = collections.Counter()
    unc = []; sources = []; reports = {}; used_names = set()
    with atomic('gene_tokens.tsv', 'gene\tchr\tkey19\tpos38\tabc\trank') as gf, \
         atomic('dist_tokens.tsv', 'lead\tkey19\tpos38') as df, atomic('lead_snps.tsv') as lf:
        for ch in CHS:
            guard(); gr, dr, lr, u, src, report, rs, co = process(ch, leads[ch], px[ch], cs[ch], genes, el, cds[ch], tissue, used_names)
            write_rows(gf, gr); write_rows(df, dr); write_rows(lf, lr)
            uk = {r[2] for r in gr}; ud = {r[1] for r in dr}; extract = sorted(uk | ud, key=lambda k: (int(k.split(':')[1]), k))
            with atomic(f'extract_chr{ch}.keys') as f: write_rows(f, ((k,) for k in extract))
            totals.update(report); totals.update(tokens_with_dup=len(gr), unique_ann=len(uk), unique_dist=len(ud), overlap=len(uk & ud),
                                                 unique_to_extract=len(extract), positive_score_tokens=sum(r[4] > 0 for r in gr))
            gc.update(r[0] for r in gr); unc.extend(u); sources.extend(src); reports[ch] = dict(report); reasons.update(rs); combos.update(co)
        need(manifest(inputs) == initial and sorted(OT.glob('*.parquet')) == files, 'Inputs changed during selection.')
        summary = dict(totals, rule='C+LDproxy0.5+SuSiE_CS; seed OT ENCODE-rE2G', cs=cssum, n_sig_nonMHC=nsig,
                       genes=len(gc), seed_tissues=SEEDS, seed_tissue_counts=tissue, element_audit=ea, cache_audit=cache_audit,
                       link_reason_counts=reasons, link_reason_combinations=combos, chromosomes=reports, input_manifest=initial,
                       missing_seed_tissues=[t for t in SEEDS if not tissue[t]['raw_rows']], cap=CAP, rank_base=0,
                       leads_with_proxy=sum(a for a, b, c in sources), leads_with_cs=sum(b for a, b, c in sources),
                       link_variants_per_lead_median=quantiles([c for a, b, c in sources], (50,))[0],
                       tokens_per_gene_q=quantiles(gc.values(), (0, 25, 50, 75, 100)), genes_at_cap=sum(n >= CAP for n in gc.values()),
                       uncapped_tokens_per_gene_q=quantiles(unc, (0, 25, 50, 75, 90, 100)),
                       frac_tokens_abc_gt0=round(totals['positive_score_tokens'] / max(1, totals['tokens_with_dup']), 4),
                       coordinates=dict(key19='hg19 full allele key', pos38='hg38 1-based',
                                        element_membership='start < pos38 <= end', gene_body='GTF inclusive +/-3000'),
                       seed_element_cache=dict(file='seed_elements.tsv.gz', scope='usable protein-coding gene elements',
                                               coordinate_system='hg38 BED', rows=ea['usable_seed_gene_elements']),
                       resource_limits=dict(virtual_memory_gib=cap, cpu_affinity_max=4, arrow_threads=1, chromosome_workers=1),
                       peak_self_rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**2,
                       peak_child_rss_gib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024**2)
        summary['status'] = 'failed_dist_shortage' if totals['dist_shortage_leads'] else ('ready' if gc else 'failed_zero_tokens')
        with atomic('selection_summary.json') as f: json.dump(summary, f, indent=2, sort_keys=True, allow_nan=False); f.write('\n')
        need(not totals['dist_shortage_leads'], 'DIST shortage: see aggregate summary; no completion marker.')
        need(gc, 'No gene tokens; no completion marker.')
    sync_directory()
    with atomic('select_seed.done') as f: json.dump(dict(status='complete', schema_version=2, gene_token_rows=totals['tokens_with_dup']), f); f.write('\n')
    sync_directory(); print(json.dumps(dict(status='complete', genes=len(gc), tokens_with_dup=totals['tokens_with_dup'])))
if __name__ == '__main__':
    try: main()
    except Exception as exc:
        print('Selection stopped: ' + (str(exc) if isinstance(exc, SelectionError) else type(exc).__name__), file=sys.stderr); sys.exit(1)
