#!/usr/bin/env python
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import sys, os, gzip, json, time, bisect, zipfile, io, resource, collections, random
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np

ROOT = _config_path("${PROJECT_ROOT}/work"); W = ROOT + "/phi_gate"
OUT = W + "/out"; PRIV = W + "/private"
UNI = ROOT + "/fset/out/uni_chr{}.tsv.gz"
VCF = ROOT + "/ref/orig_index/chr{}.vcf.gz"
GTF = ROOT + "/ref/deductive/gencode.nochr.gtf.gz"
CCRE = ROOT + "/ref/b6_cards/ccre.s.bed"; RMSK = ROOT + "/ref/repeats/rmsk.txt.gz"; CPG = ROOT + "/ref/cpg/cpgIslandExt.txt.gz"
BW = ROOT + "/ref/mappability/k36.Umap.MultiTrackMappability.bw"
REMAP_BY_CHR = ROOT + "/fset/out/remap_by_chr/chr{}.bed"; TOPTF = ROOT + "/fset/out/remap_top_tfs.txt"
GPN = ROOT + "/ref/gpnmsa/scores.tsv.bgz"
BBJ_ZIP = ROOT + "/ref/bbj/hum0197.v3.BBJ.Hei.v1.zip"
BBJ_MEMBER = "hum0197.v3.BBJ.Hei.v1/GWASsummary_Height_Japanese_SakaueKanai2020.auto.txt.gz"
LEADS = ROOT + "/prs/out/hei/lead_snps.tsv"
R2_MIN = 0.3; FLANK = 3000; TEAM_R2 = 0.8
SEED = 20260928; NSAMP = 10000
BLOCK = 2000
CHUNK_BP = 2_000_000; WORKERS = 8
KNOTS = (0.2, 0.4, 0.6, 0.8)
MIN_TEAMS = 3
NPERM = 1000; VAR_SHARE = 0.5; Q95 = 95; FDR = 0.10
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

def log(stage, **kw):
    kw.update(stage=stage, t=time.strftime("%Y-%m-%dT%H:%M:%S"), rss_gb=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6, 2))
    print(json.dumps(kw), flush=True)

def atomic_json(obj, path):
    tmp = path + ".tmp"; json.dump(obj, open(tmp, "w"), indent=1); os.replace(tmp, path)

def merge(iv):
    iv = sorted(iv); m = []
    for s, e, g in iv:
        m.append((s, e, g))
    return m

def load_universe(N):
    keys, p38, maf = [], [], []
    with gzip.open(UNI.format(N), "rt") as f:
        h = f.readline().rstrip("\n").split("\t"); ix = {c: i for i, c in enumerate(h)}
        for l in f:
            a = l.rstrip("\n").split("\t")
            if a[ix["is_coding"]] != "0" or float(a[ix["r2"]]) < R2_MIN: continue
            k = a[ix["key"]]
            if "," in k.split(":")[3]: raise SystemExit("multi-ALT key in universe")
            keys.append(k); p38.append(int(a[ix["pos38"]])); maf.append(float(a[ix["maf"]]))
    assert len(set(keys)) == len(keys), "duplicate universe keys"
    return keys, np.array(p38), np.array(maf)

def load_genes(N):
    g = []
    for l in gzip.open(GTF, "rt"):
        if l[0] == "#": continue
        a = l.split("\t")
        if a[2] != "gene" or 'gene_type "protein_coding"' not in a[8]: continue
        if a[0].replace("chr", "") != str(N): continue
        gid = a[8].split('gene_id "')[1].split('"')[0].split(".")[0]
        g.append((max(1, int(a[3]) - FLANK), int(a[4]) + FLANK, gid))
    g.sort(key=lambda x: x[2]); return g

def assign(p38, genes):
    order = np.argsort(p38); ps = p38[order]
    members = []
    for s, e, gid in genes:
        i, j = np.searchsorted(ps, s, "left"), np.searchsorted(ps, e, "right")
        members.append(np.sort(order[i:j]))
    return members

def _team_block(D):
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    m = D.shape[1]
    if m == 1: return np.zeros(1, int)
    mu = np.nanmean(D, 0); D = np.where(np.isnan(D), mu, D); D = D - D.mean(0)
    sd = np.sqrt((D ** 2).sum(0)); poly = sd > 0
    Z = np.zeros_like(D); Z[:, poly] = D[:, poly] / sd[poly]
    R2 = (Z.T @ Z) ** 2; np.fill_diagonal(R2, 1.0)
    R2[~poly, :] = 0; R2[:, ~poly] = 0; R2[~poly, ~poly] = 0; np.fill_diagonal(R2, 1.0)
    dist = np.clip(1 - R2, 0, 1); np.fill_diagonal(dist, 0)
    L = linkage(squareform(dist, checks=False), method="complete")
    lab = fcluster(L, np.nextafter(1 - TEAM_R2, -np.inf), criterion="distance") - 1
    for t in np.unique(lab):
        idx = np.where(lab == t)[0]
        if len(idx) > 1: assert (R2[np.ix_(idx, idx)] > TEAM_R2).all(), "team violates r2 rule"
    return lab

def _worker(args):
    N, lo, hi, want, samples = args
    from cyvcf2 import VCF as V
    v = V(VCF.format(N), samples=samples, lazy=True, threads=1)
    ctg = [c for c in v.seqnames if c.replace("chr", "") == str(N)]; assert len(ctg) == 1, "contig"
    got_k, got_d = [], []
    for rec in v("%s:%d-%d" % (ctg[0], lo, hi)):
        for ai, alt in enumerate(rec.ALT):
            k = "chr%s:%d:%s:%s" % (N, rec.POS, rec.REF, alt)
            if k not in want: continue
            ds = rec.format("DS")
            if ds is None: raise RuntimeError("DS missing")
            col = ds[:, ai] if ds.shape[1] > ai else ds[:, 0]
            got_k.append(k); got_d.append(col.astype(np.float32))
    return got_k, got_d

def build_teams(N, keys, members_all):
    from cyvcf2 import VCF as V
    names = V(VCF.format(N)).samples
    rng = random.Random(SEED); samples = sorted(rng.sample(range(len(names)), NSAMP)); samples = [names[i] for i in samples]
    idx = sorted(set(np.concatenate(members_all).tolist()), key=lambda i: int(keys[i].split(":")[1]))
    pos19 = [int(keys[i].split(":")[1]) for i in idx]; want = {keys[i]: i for i in idx}
    jobs = []; lo = pos19[0]
    while lo <= pos19[-1]:
        hi = lo + CHUNK_BP - 1
        sub = {k: i for k, i in want.items() if lo <= int(k.split(":")[1]) <= hi}
        if sub: jobs.append((N, lo, hi, sub, samples))
        lo = hi + 1
    from multiprocessing import Pool
    team = np.full(len(keys), -1, int); nxt = 0; done = 0
    with Pool(WORKERS) as pl:
        for gk, gd in pl.imap(_worker, jobs):
            done += 1
            order = np.argsort([int(k.split(":")[1]) for k in gk], kind="stable")
            gk = [gk[o] for o in order]; gd = [gd[o] for o in order]
            for b in range(0, len(gk), BLOCK):
                D = np.stack(gd[b:b + BLOCK], 1)
                lab = _team_block(D)
                for k, t in zip(gk[b:b + BLOCK], lab): team[want[k]] = nxt + t
                nxt += lab.max() + 1
            if done % 5 == 0 or done == len(jobs): log("teams", chunks_done=done, chunks_total=len(jobs), teams=int(nxt))
    miss = sum(1 for i in idx if team[i] < 0)
    if miss: raise SystemExit("dosage not found for %d assigned variants" % miss)
    return team

def overlap_dist(iv, P):
    starts = [a for a, _, _ in iv]; maxend = np.maximum.accumulate([b for _, b, _ in iv]) if iv else []
    out = []
    for p in P:
        p0 = p - 1; i = bisect.bisect_right(starts, p0) - 1; val = None; j = i
        while j >= 0 and maxend[j] > p0:
            a, b, vv = iv[j]
            if a <= p0 < b: val = vv; break
            j -= 1
        if val is None:
            c = []
            if i >= 0: c.append(p0 - iv[i][1] + 1)
            k = i + 1
            if k < len(iv): c.append(iv[k][0] - p0)
            out.append((None, min(c) if c else None))
        else: out.append((val, 0))
    return out

def load_bed(path, cc, sc, ec, vc, gz, chrom):
    iv = []
    with (gzip.open(path, "rt") if gz else open(path)) as f:
        for l in f:
            a = l.rstrip("\n").split("\t")
            if len(a) <= max(cc, sc, ec, vc or 0) or a[cc] != chrom: continue
            iv.append((int(a[sc]), int(a[ec]), a[vc] if vc is not None else "1"))
    iv.sort(); return iv

def annotations(N, keys, p38, sel):
    chrom = "chr%s" % N; P = [int(p38[i]) for i in sel]; cols = {}
    for name, (path, cc, sc, ec, vc, gz) in {"ccre": (CCRE, 0, 1, 2, 3, False), "rep": (RMSK, 5, 6, 7, 11, True), "cpg": (CPG, 1, 2, 3, None, True)}.items():
        r = overlap_dist(load_bed(path, cc, sc, ec, vc, gz, chrom), P)
        if name == "cpg": cols["cpg"] = np.array([1.0 if v else 0.0 for v, _ in r])
        else:
            cls = [v or "none" for v, _ in r]
            for c in sorted(set(cls)): cols["%s_class=%s" % (name, c)] = np.array([1.0 if x == c else 0.0 for x in cls])
        cols[name + "_dist"] = np.array([np.nan if d is None else float(d) for _, d in r])
    import pyBigWig
    bw = pyBigWig.open(BW); cols["map_k36"] = np.array([np.nan if (x is None or x != x) else float(x) for x in (bw.values(chrom, p - 1, p)[0] for p in P)]); bw.close()
    top = [l.strip() for l in open(TOPTF)]
    ivs = sorted((int(a[1]), int(a[2]), a[3].split(":")[0]) for a in (l.split("\t", 5) for l in open(REMAP_BY_CHR.format(N))) if a[0] == chrom)
    st = [a for a, _, _ in ivs]; hits = []
    for p in P:
        p0 = p - 1; j = bisect.bisect_right(st, p0) - 1; s = set()
        while j >= 0 and p0 - st[j] < 20000:
            a, b, tf = ivs[j]
            if a <= p0 < b: s.add(tf)
            j -= 1
        hits.append(s)
    cols["tf_n"] = np.array([float(len(s)) for s in hits])
    for tf in top: cols["tf_" + tf] = np.array([1.0 if tf in s else 0.0 for s in hits])
    import pysam
    tb = pysam.TabixFile(GPN); ctg = [c for c in tb.contigs if c.replace("chr", "") == str(N)][0]
    g = np.full(len(sel), np.nan)
    for n, i in enumerate(sel):
        _, _, ref, alt = keys[i].split(":")
        if len(ref) != 1 or len(alt) != 1: continue
        for row in tb.fetch(ctg, int(p38[i]) - 1, int(p38[i])):
            a = row.split("\t"); r2, a2 = a[2].upper(), a[3].upper()
            if (r2, a2) in ((ref, alt), (COMP[ref], COMP[alt])): g[n] = float(a[4])
            elif (r2, a2) in ((alt, ref), (COMP[alt], COMP[ref])): g[n] = -float(a[4])
    cols["gpn_msa"] = g
    return cols

def bbj_absbeta(N, keys):
    kidx = {}
    for i, k in enumerate(keys):
        _, p, r, a = k.split(":"); kidx.setdefault(int(p), []).append((r, a, i))
    out = np.full(len(keys), np.nan); stats = collections.Counter()
    zf = zipfile.ZipFile(BBJ_ZIP)
    f = io.TextIOWrapper(gzip.GzipFile(fileobj=zf.open(BBJ_MEMBER)))
    h = f.readline().rstrip("\n").split("\t"); ix = {c: j for j, c in enumerate(h)}
    for l in f:
        a = l.rstrip("\n").split("\t")
        if a[ix["CHR"]] != str(N): continue
        cand = kidx.get(int(a[ix["BP"]]))
        if not cand: continue
        a1, a0 = a[ix["ALLELE1"]], a[ix["ALLELE0"]]
        for r, alt, i in cand:
            if len(r) == 1 and len(alt) == 1 and {r, alt} in ({"A", "T"}, {"C", "G"}): stats["palindromic_skipped"] += 1; continue
            ok = {a1, a0} == {r, alt} or (len(r) == 1 and len(alt) == 1 and {COMP.get(a1), COMP.get(a0)} == {r, alt})
            if not ok: continue
            if not np.isnan(out[i]): raise SystemExit("multiple BBJ rows for one variant")
            out[i] = abs(float(a[ix["BETA"]])); stats["matched"] += 1
    return out, dict(stats)

def lead_dist(N, p38):
    L = np.array(sorted(int(a[4]) for a in (l.rstrip("\n").split("\t") for l in open(LEADS)) if a[0] == str(N)))
    if len(L) == 0: return np.full(len(p38), np.nan)
    j = np.clip(np.searchsorted(L, p38), 1, len(L) - 1)
    return np.minimum(np.abs(p38 - L[j - 1]), np.abs(L[j] - p38)).astype(float)

def run_chr(N):
    od = "%s/chr%s" % (OUT, N); pdir = "%s/chr%s" % (PRIV, N)
    if os.path.exists(od + "/chr.done"): raise SystemExit("done marker exists: " + od)
    os.makedirs(od, exist_ok=True); os.makedirs(pdir, mode=0o700, exist_ok=True)
    keys, p38, maf = load_universe(N); log("universe", chr=N, n=len(keys))
    genes = load_genes(N); members = assign(p38, genes)
    keep = [g for g in range(len(genes)) if len(members[g]) > 0]
    log("genes", chr=N, genes=len(genes), genes_with_variants=len(keep), assignments=int(sum(len(m) for m in members)))
    team = build_teams(N, keys, [members[g] for g in keep])
    sel = np.array(sorted(set(np.concatenate([members[g] for g in keep]).tolist())))
    ann = annotations(N, keys, p38, sel); log("annotations", chr=N, n_cols=len(ann))
    absb, bst = bbj_absbeta(N, keys); ld = lead_dist(N, p38); log("bbj", chr=N, **bst)
    pos_in_sel = {int(v): n for n, v in enumerate(sel)}
    inst_gene, inst_team, inst_members = [], [], []
    for g in keep:
        by = collections.defaultdict(list)
        for i in members[g]: by[int(team[i])].append(pos_in_sel[int(i)])
        for t, mem in sorted(by.items()):
            inst_gene.append(g); inst_team.append(t); inst_members.append(np.array(mem))
    names = sorted(ann); A = np.column_stack([ann[c] for c in names])
    X = np.vstack([np.nanmean(A[m], 0) if np.isfinite(A[m]).any() else np.full(len(names), np.nan) for m in inst_members]) if inst_members else np.zeros((0, len(names)))
    lens = np.array([len(m) for m in inst_members]); flat = np.concatenate(inst_members)
    np.savez(pdir + "/chr_build.tmp.npz", X=X, gene=np.array(inst_gene), team=np.array(inst_team), mem_len=lens, mem_flat=flat,
             maf=maf[sel], absb=absb[sel], lead=ld[sel], names=np.array(names))
    os.replace(pdir + "/chr_build.tmp.npz", pdir + "/chr_build.npz")
    diag = dict(chr=N, universe=len(keys), variants_in_genes=int(len(sel)), genes=len(keep), team_instances=len(inst_gene),
                teams_unique=int(len(set(team[sel].tolist()))), bbj=bst, bbj_match_rate=round(bst.get("matched", 0) / max(1, len(sel)), 4),
                annot_coverage={c: round(float(np.isfinite(ann[c]).mean()), 4) for c in names})
    atomic_json(diag, od + "/diagnostics.json"); open(od + "/chr.done", "w").write("ok %s\n" % time.strftime("%Y-%m-%dT%H:%M:%S"))
    log("chr_done", chr=N)

def ns_basis(x, knots):
    k = np.asarray(knots, float); K = len(k)
    d = lambda j: (np.clip(x - k[j], 0, None) ** 3 - np.clip(x - k[-1], 0, None) ** 3) / (k[-1] - k[j])
    cols = [x] + [d(j) - d(K - 2) for j in range(K - 2)]
    return np.column_stack(cols)

def rank_avg(v):
    from scipy.stats import rankdata
    return rankdata(v)

def run_pool(chrs):
    od = OUT + "/pooled_" + "_".join(map(str, chrs))
    if os.path.exists(od + "/gate.done"): raise SystemExit("done marker exists")
    os.makedirs(od, exist_ok=True)
    Bs = [np.load("%s/chr%s/chr_build.npz" % (PRIV, N), allow_pickle=False) for N in chrs]
    names = list(Bs[0]["names"]); assert all(list(b["names"]) == names for b in Bs[1:]), "annotation columns differ across chromosomes"
    maf = np.concatenate([b["maf"] for b in Bs]); ab = np.concatenate([b["absb"] for b in Bs]); ld = np.concatenate([b["lead"] for b in Bs])
    ok = np.isfinite(ab) & (ab > 0) & (maf > 0) & np.isfinite(ld)
    x1 = np.log10(maf); x2 = np.log10(ld + 1)
    def kn(x): return [x[ok].min()] + [np.quantile(x[ok], q) for q in KNOTS] + [x[ok].max()]
    Xs = np.column_stack([np.ones(len(maf)), ns_basis(x1, kn(x1)), ns_basis(x2, kn(x2))])
    beta, *_ = np.linalg.lstsq(Xs[ok], np.log(ab[ok]), rcond=None)
    res = np.full(len(maf), np.nan); res[ok] = np.log(ab[ok]) - Xs[ok] @ beta
    log("residual_fit", n_fit=int(ok.sum()), n_total=len(maf), rank=int(np.linalg.matrix_rank(Xs[ok])), n_params=Xs.shape[1])
    Xall, G, Y = [], [], []; off = 0; goff = 0
    for b in Bs:
        L = b["mem_len"]; F = b["mem_flat"] + off; st = np.concatenate([[0], np.cumsum(L)])
        for n in range(len(L)):
            r = res[F[st[n]:st[n + 1]]]; Y.append(np.nanmax(r) if np.isfinite(r).any() else np.nan)
        Xall.append(b["X"]); G.append(b["gene"] + goff); off += len(b["maf"]); goff += int(b["gene"].max()) + 1
    X = np.vstack(Xall); G = np.concatenate(G); Y = np.array(Y)
    m1 = {}
    for c, name in enumerate(names):
        x = X[:, c]; f = np.isfinite(x)
        if f.sum() < 2: m1[name] = None; continue
        xs, gs = x[f], G[f]; tot = ((xs - xs.mean()) ** 2).sum()
        if tot == 0: m1[name] = 0.0; continue
        gm = np.bincount(gs, weights=xs) / np.maximum(np.bincount(gs), 1)
        m1[name] = float(((xs - gm[gs]) ** 2).sum() / tot)
    cand = [n for n in names if m1[n] is not None and m1[n] >= VAR_SHARE]
    log("measure1", n_annot=len(names), n_candidates=len(cand))
    genes = np.unique(G); rows = []
    for c_i, name in enumerate(cand):
        c = names.index(name); obs = []; null = np.zeros(NPERM); ng = 0
        rng = np.random.default_rng([SEED, c])
        for g in genes:
            ii = np.where((G == g) & np.isfinite(X[:, c]) & np.isfinite(Y))[0]
            if len(ii) < MIN_TEAMS: continue
            rx = rank_avg(X[ii, c]); ry = rank_avg(Y[ii])
            if rx.std() == 0 or ry.std() == 0: continue
            zx = (rx - rx.mean()) / rx.std(); zy = (ry - ry.mean()) / ry.std()
            obs.append(float((zx * zy).mean()))
            P = np.argsort(rng.random((NPERM, len(ii))), 1)
            null += (zy[P] * zx).mean(1); ng += 1
        if ng == 0: rows.append(dict(annot=name, genes=0, status="not_estimable", p=1.0)); continue
        o = float(np.mean(obs)); nm = null / ng
        rows.append(dict(annot=name, genes=ng, obs=o, null_q95=float(np.percentile(nm, Q95)), p=float((1 + (nm >= o).sum()) / (NPERM + 1)), status="ok"))
        if (c_i + 1) % 10 == 0: log("measure2", done=c_i + 1, total=len(cand))
    p = np.array([r["p"] for r in rows]); m = len(p)
    if m:
        o = np.argsort(p); q = np.empty(m); q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
        for r, qq in zip(rows, q):
            r["q"] = float(min(qq, 1.0)); r["pass"] = bool(r["status"] == "ok" and r["obs"] > r["null_q95"] and r["q"] < FDR)
    passed = [r["annot"] for r in rows if r.get("pass")]
    verdict = "PROCEED_TO_PHI_TRAINING" if passed else "NO_DISCRIMINABLE_EFFECT_HEIGHT_CHR" + "+".join(map(str, chrs))
    atomic_json(dict(chrs=chrs, measure1=m1, candidates=cand, measure2=rows, passed=passed, verdict=verdict,
                     instances=int(len(G)), genes=int(len(genes)), residual_fit_n=int(ok.sum())), od + "/decision.json")
    open(od + "/gate.done", "w").write(verdict + "\n"); log("done", verdict=verdict, passed=len(passed))

if __name__ == "__main__":
    os.nice(10)
    if sys.argv[1] == "chr": run_chr(sys.argv[2])
    elif sys.argv[1] == "pool": run_pool(sys.argv[2:])
    else: raise SystemExit(__doc__)
