#!/usr/bin/env python
import os as _os
if not _os.environ.get("PROJECT_ROOT", "").strip():
    raise ValueError("PROJECT_ROOT must be set and nonblank")
import gzip, json, os, sys, bisect, collections
CHR = sys.argv[1]; SAMPLE = sys.argv[2]
ROOT = (_os.environ["PROJECT_ROOT"] + '/work'); OUT = f"{ROOT}/fset/out"
CCRE = f"{ROOT}/ref/b6_cards/ccre.s.bed"
RMSK = f"{ROOT}/ref/repeats/rmsk.txt.gz"
CPG = f"{ROOT}/ref/cpg/cpgIslandExt.txt.gz"
BW = f"{ROOT}/ref/mappability/k36.Umap.MultiTrackMappability.bw"
REMAP = f"{ROOT}/ref/remap/remap2022_nr_macs2_hg38_v1_0.bed.gz"
RE2G_DIR = f"{ROOT}/ref/re2g"
TOPTF = f"{OUT}/remap_top_tfs.txt"
chrname = f"chr{CHR}"

want = {}
with open(SAMPLE) as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        if p[0] == "key" or len(p) < 4: continue
        if p[1] != CHR: continue
        want[p[0]] = (int(p[2]), p[3])
if not want:
    print(json.dumps({"chr": CHR, "variants": 0})); sys.exit(0)
keys = sorted(want, key=lambda k: want[k][0])
P = [want[k][0] for k in keys]

def load_bed(path, chrom_col, start_col, end_col, val_col, gz=True, want_chrom=chrname, skip_prefix=None):
    iv = []
    op = gzip.open(path, "rt") if gz else open(path)
    with op as f:
        for line in f:
            if skip_prefix and line.startswith(skip_prefix): continue
            p = line.rstrip("\n").split("\t")
            if len(p) <= max(chrom_col, start_col, end_col, val_col if val_col is not None else 0): continue
            if p[chrom_col] != want_chrom: continue
            iv.append((int(p[start_col]), int(p[end_col]), p[val_col] if val_col is not None else "1"))
    iv.sort()
    return iv

def overlap_and_dist(iv, positions):
    if not iv: return [(None, None)] * len(positions)
    starts = [a for a, _, _ in iv]; ends = [b for _, b, _ in iv]
    maxend = []
    m = -1
    for b in ends:
        m = max(m, b); maxend.append(m)
    out = []
    for p in positions:
        p0 = p - 1
        i = bisect.bisect_right(starts, p0) - 1
        val = None; dist = None
        j = i
        while j >= 0 and maxend[j] > p0:
            a, b, v = iv[j]
            if a <= p0 < b: val = v; break
            j -= 1
        if val is None:
            cands = []
            if i >= 0: cands.append(p0 - iv[i][1] + 1)
            k = bisect.bisect_right(starts, p0)
            if k < len(iv): cands.append(iv[k][0] - p0)
            dist = min(c for c in cands) if cands else None
        else:
            dist = 0
        out.append((val, dist))
    return out

res = {k: {} for k in keys}
cc = overlap_and_dist(load_bed(CCRE, 0, 1, 2, 3, gz=False), P)
for k, (v, d) in zip(keys, cc): res[k]["ccre_class"] = v or "none"; res[k]["ccre_dist"] = d
rm = overlap_and_dist(load_bed(RMSK, 5, 6, 7, 11), P)
for k, (v, d) in zip(keys, rm): res[k]["rep_class"] = v or "none"; res[k]["rep_dist"] = d
cg = overlap_and_dist(load_bed(CPG, 1, 2, 3, None), P)
for k, (v, d) in zip(keys, cg): res[k]["cpg"] = 1 if v else 0; res[k]["cpg_dist"] = d
try:
    import pyBigWig
    bwf = pyBigWig.open(BW)
    for k in keys:
        p = want[k][0]
        try:
            v = bwf.values(chrname, p - 1, p)[0]
        except Exception:
            v = None
        res[k]["map_k36"] = None if v is None or v != v else float(v)
    bwf.close()
    map_status = "ok"
except Exception as e:
    for k in keys: res[k]["map_k36"] = None
    map_status = f"unavailable:{type(e).__name__}"
top = [l.strip() for l in open(TOPTF)] if os.path.exists(TOPTF) else []
topset = set(top)
tfhit = collections.defaultdict(set)
ivs = []
percr = f"{ROOT}/fset/out/remap_by_chr/{chrname}.bed"
src = percr if os.path.exists(percr) else REMAP
op = open(src) if src == percr else gzip.open(REMAP, "rt")
with op as f:
    for line in f:
        p = line.split("\t", 5)
        if p[0] != chrname: continue
        ivs.append((int(p[1]), int(p[2]), p[3].split(":")[0]))
ivs.sort()
starts = [a for a, _, _ in ivs]
for k in keys:
    p0 = want[k][0] - 1
    i = bisect.bisect_right(starts, p0)
    j = i - 1
    seen = set()
    while j >= 0 and p0 - starts[j] < 20000:
        a, b, tf = ivs[j]
        if a <= p0 < b: seen.add(tf)
        j -= 1
    tfhit[k] = seen
for k in keys:
    s = tfhit[k]
    res[k]["tf_n"] = len(s)
    for tf in top: res[k][f"tf_{tf}"] = 1 if tf in s else 0
re2g_files = []
cols = ["key", "chr", "pos38", "domain"] + sorted({c for k in keys for c in res[k]})
tmp = f"{_os.environ['PROJECT_ROOT']}/work/prs/out/hei_seed/feat/feat_chr{CHR}.tsv.gz.tmp"
with gzip.open(tmp, "wt") as out:
    out.write("\t".join(cols) + "\n")
    for k in keys:
        row = {"key": k, "chr": CHR, "pos38": want[k][0], "domain": want[k][1]}
        row.update(res[k])
        out.write("\t".join("" if row.get(c) is None else str(row.get(c, "")) for c in cols) + "\n")
os.replace(tmp, f"{_os.environ['PROJECT_ROOT']}/work/prs/out/hei_seed/feat/feat_chr{CHR}.tsv.gz")
print(json.dumps({"chr": CHR, "variants": len(keys), "n_cols": len(cols), "mappability": map_status,
                  "re2g_biosamples": len(re2g_files), "top_tf_used": len(top),
                  "ccre_overlap": sum(1 for k in keys if res[k]["ccre_class"] != "none"),
                  "tf_any": sum(1 for k in keys if res[k]["tf_n"] > 0)}))
