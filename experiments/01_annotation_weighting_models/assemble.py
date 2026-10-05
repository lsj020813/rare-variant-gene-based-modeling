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


import sys, os, json
def require(cond, msg):
    if not cond: raise RuntimeError(f"{msg}")
require(len(sys.argv) == 2 and sys.argv[1].isdigit() and 1 <= int(sys.argv[1]) <= 22 and sys.argv[1] == str(int(sys.argv[1])), "usage: assemble.py <1..22>")
N = sys.argv[1]
R = _config_path("${PROJECT_ROOT}/work/ref15")
GP  = f"{R}/groupfiles_bwg/chr{N}.B_3kb_re2g.txt"
FAV = f"{R}/annot/extract/chr{N}.annot.tsv"
GPN = f"{R}/annot/gpn/chr{N}.gpn.tsv"
T2  = f"{R}/annot/t2/chr{N}.t2.tsv"
OUT = f"{R}/annot/fm"; os.makedirs(OUT, exist_ok=True)
O = f"{OUT}/chr{N}.fm.tsv"
for p in (GP, FAV, GPN, T2):
    require(os.path.getsize(p) > 0, f"missing {p}")
import gzip, re as _re
GTF = f"{R}/deductive/gencode.sorted.gtf.gz"
sym = {}
with gzip.open(GTF, "rt") as fh:
    for line in fh:
        if line[0] == "#": continue
        f = line.split("\t")
        if f[2] != "gene" or f[0] != f"chr{N}": continue
        i = f[8].find('gene_id "'); gid = f[8][i+9:f[8].find('"', i+9)].split(".")[0]
        j = f[8].find('gene_name "'); sym[gid] = f[8][j+11:f[8].find('"', j+11)] if j >= 0 else ""
require(sym, "no genes in GTF for this chrom")
_num = _re.compile(r"^-?[0-9.]+([eE][-+]?[0-9]+)?$")
def parse_gh(txt, gene_sym):
    if not txt: return ("", "", "")
    m = _re.search(r"Name=([0-9.]+)", txt); elem = m.group(1) if m else ""
    links = _re.findall(r"connected_gene=([^;]+);score=([0-9.]+)", txt)
    link = next((sc for g, sc in links if g == gene_sym), "")
    return (elem, str(len(links)), link)

EXPECTED = {"FAV": ["key37","cons","epi_active","epi_repr","epi_trans","tf","cage_prom","genehancer","linsight"],
            "GPN": ["key37","gpn_msa"]}
def load(path, which):
    d = {}
    with open(path, encoding="utf-8", newline="") as f:
        hdr = next(f).rstrip("\n").split("\t")
        require(hdr == EXPECTED[which], f"{which} header mismatch: {hdr}")
        for line in f:
            a = line.rstrip("\n").split("\t")
            require(len(a) == len(hdr), f"{path} row width {len(a)} != header {len(hdr)}")
            require(a[0] not in d, f"duplicate key in {path}")
            d[a[0]] = a[1:]
    return hdr[1:], d

fav_cols, fav = load(FAV, "FAV")
gpn_cols, gpn = load(GPN, "GPN")
t2 = {}
with open(T2) as f:
    t2_cols = next(f).rstrip("\n").split("\t")[2:]
    for line in f:
        a = line.rstrip("\n").split("\t")
        require(len(a) == len(t2_cols) + 2, "t2 row width")
        t2[(a[0], a[1])] = a[2:]
require(fav_cols[5] == "cage_prom" and fav_cols[6] == "genehancer", f"fav cols {fav_cols}")
num_cols = [c for c in fav_cols if c not in ("cage_prom", "genehancer")] + gpn_cols
cols = num_cols
derived = ["is_cage_prom", "gh_elem_score", "gh_n_genes", "gh_link_score"]
vlen = {}
with open(GP) as f:
    for line in f:
        a = line.split()
        if len(a) < 3: continue
        if a[1] == "var": vlen[a[0]] = len(a)
        elif a[1] == "anno": require(vlen.get(a[0]) == len(a), f"var/anno cardinality mismatch for a gene")
        else: require(False, f"unknown record tag {a[1]!r}")
pairs = 0; dup = set(); nas = {c: 0 for c in cols}
OT = O + ".tmp"
for p in (O + ".done",): 
    if os.path.exists(p): os.remove(p)
with open(OT, "w", encoding="utf-8", newline="") as o:
    o.write("key37\tgene\tgene_base\t" + "\t".join(cols) + "\t" + "\t".join(c + "_na" for c in cols) + "\t" + "\t".join(derived) + "\t" + "\t".join(t2_cols) + "\n")
    with open(GP) as f:
        for line in f:
            a = line.split()
            if len(a) < 3 or a[1] != "var": continue
            gene = a[0]
            for v in a[2:]:
                k = (v, gene)
                require(k not in dup, f"dup pair {k}")
                dup.add(k); pairs += 1
                fv = fav.get(v); gv = gpn.get(v)
                fvd = dict(zip(fav_cols, fv)) if fv else {}
                vals = [fvd.get(c, "") for c in fav_cols if c not in ("cage_prom", "genehancer")] + (gv if gv else [""] * len(gpn_cols))
                for x in vals: require(x == "" or _num.match(x), f"non-numeric value in {v}")
                cage = "1" if fvd.get("cage_prom") else "0"
                gh = parse_gh(fvd.get("genehancer", ""), sym.get(gene.split(".")[0], ""))
                der = [cage, *gh]
                na = ["1" if x == "" else "0" for x in vals]
                for c, x in zip(cols, vals):
                    if x == "": nas[c] += 1
                tv = t2.get(k); require(tv is not None, f"pair {k} missing in T2")
                o.write(f"{v}\t{gene}\t{gene.split('.')[0]}\t" + "\t".join(vals) + "\t" + "\t".join(na) + "\t" + "\t".join(der) + "\t" + "\t".join(tv) + "\n")
require(pairs > 0, "0 pairs")
SUM = f"{R}/groupfiles_bwg/chr{N}.summary.json"
if os.path.exists(SUM):
    sm = json.load(open(SUM))
    require(sm["pairs"] == pairs, f"pairs {pairs} != summary {sm['pairs']}")
    require(sm["distinct"] == len({v for v, _ in dup}), "distinct variants != summary")
stats = {"chr": N, "pairs": pairs, "derived": derived, "genes": len({g for _, g in dup}),
         "variants": len({v for v, _ in dup}),
         "na_pct": {c: round(nas[c] / pairs * 100, 1) for c in cols}}
vset = {v for v, _ in dup}
for name, src, mn in (("FAVOR", fav, 0.85), ("GPN", gpn, 0.85)):
    rate = len(vset & set(src)) / len(vset)
    stats[f"hit_{name}"] = round(rate, 4)
    require(rate >= mn, f"{name} key-match rate {rate:.3f} < {mn}")
with open(f"{O}.stats.json.tmp", "w") as fh: json.dump(stats, fh)
os.replace(OT, O); os.replace(f"{O}.stats.json.tmp", f"{O}.stats.json")
with open(f"{O}.done", "w") as fh: fh.write(f"ok {pairs}\n")
print("FMSTATS", json.dumps(stats))
print(f"[chr{N}] ASSEMBLE_DONE")
