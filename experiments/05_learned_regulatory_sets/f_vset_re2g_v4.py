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


import gzip, os, sys, json, bisect, collections
ROOT = _config_path("${PROJECT_ROOT}/work"); OUT = f"{ROOT}/fset/primary"; RE2G_DIR = f"{ROOT}/ref/re2g"
S = sys.argv[1]
by_chr = collections.defaultdict(dict)
with open(S) as f:
    hdr = f.readline().rstrip("\n").split("\t"); ix = {h: i for i, h in enumerate(hdr)}
    for line in f:
        p = line.rstrip("\n").split("\t"); by_chr[p[ix["chr"]]][p[ix["key"]]] = int(p[ix["pos38"]])
files = sorted(f for f in os.listdir(RE2G_DIR) if f.endswith(".bed.gz"))
bios = {}
data = {bi: collections.defaultdict(list) for bi in range(len(files))}
for bi, fn in enumerate(files):
    with gzip.open(os.path.join(RE2G_DIR, fn), "rt") as f:
        for line in f:
            if line.startswith("#"): continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 30: continue
            try: sc = float(p[29])
            except ValueError: continue
            bios.setdefault(bi, p[27])
            data[bi][p[0]].append((int(p[1]), int(p[2]), sc, p[8], p[13], p[4], p[15]))
with open(f"{OUT}/vset_re2g_biosamples.tsv", "w") as o:
    o.write("index\tfile\tbiosample\n")
    for bi, fn in enumerate(files): o.write(f"{bi}\t{fn}\t{bios.get(bi,'')}\n")
def hits(iv, pos):
    starts = iv[0]; recs = iv[1]; p0 = pos - 1
    i = bisect.bisect_right(starts, p0) - 1; best = None
    while i >= 0 and p0 - starts[i] < 5000:
        a, b = recs[i][0], recs[i][1]
        if a <= p0 < b and (best is None or recs[i][2] > best[2]): best = recs[i]
        i -= 1
    return best
summary = {}
for ch in sorted(by_chr, key=lambda c: int(c)):
    want = by_chr[ch]; keys = sorted(want, key=lambda k: want[k]); chrname = f"chr{ch}"
    prepared = {}
    for bi in data:
        iv = sorted(data[bi].get(chrname, [])); prepared[bi] = ([a for a, *_ in iv], iv)
    cols = []
    for bi in range(len(files)): cols += [f"re2g{bi}_score", f"re2g{bi}_target", f"re2g{bi}_target_ensg", f"re2g{bi}_class", f"re2g{bi}_dist_tss"]
    tmp = f"{OUT}/vset_re2g_chr{ch}.tsv.gz.tmp"; nany = 0
    with gzip.open(tmp, "wt") as o:
        o.write("\t".join(["key", "chr", "pos38"] + cols) + "\n")
        for k in keys:
            row = [k, ch, str(want[k])]; any_ = False
            for bi in range(len(files)):
                h = hits(prepared[bi], want[k]) if prepared[bi][0] else None
                if h: row += [f"{h[2]:.6g}", h[3], h[4], h[5], h[6]]; any_ = True
                else: row += ["", "", "", "", ""]
            nany += any_; o.write("\t".join(row) + "\n")
    os.replace(tmp, f"{OUT}/vset_re2g_chr{ch}.tsv.gz")
    summary[ch] = {"n": len(keys), "any_re2g": nany}; print(json.dumps({ch: summary[ch]}), flush=True)
json.dump({"biosamples": bios, "per_chr": summary}, open(f"{OUT}/vset_re2g.meta.json", "w"), indent=1)
open(f"{OUT}/vset_re2g.done", "w").write("ok\n")
