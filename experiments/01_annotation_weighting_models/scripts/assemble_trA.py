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


import sys, os, gzip, json, re, time, argparse
import numpy as np
R = _config_path("${PROJECT_ROOT}/work/ref"); A = f"{R}/annot"; L1 = f"{A}/bbj_l1"
OUTDIR = _config_path("${PROJECT_ROOT}/work/run_trackA/annot_trA")
ap = argparse.ArgumentParser(); ap.add_argument("--chrs", required=True); ap.add_argument("--out", default="trA_annot_missing.tsv.gz")
args = ap.parse_args(); CHRS = [c for c in args.chrs.split(",") if c]
T0 = time.time()
def log(m): print(f"[asm {time.time()-T0:6.0f}s] {m}", flush=True)
def require(c, m):
    if not c: raise RuntimeError("GATE FAIL: " + m)
def done(p): require(os.path.exists(p + ".done") and os.path.getsize(p + ".done") > 0, f"not done: {p}")
CACHE_COLS = json.load(open(f"{A}/cache/fm_all.stats.json"))["cols"]; require(len(CACHE_COLS) == 33, "cache cols")
COLS = CACHE_COLS + ["cadd"]
T1 = ["cons", "epi_active", "epi_repr", "epi_trans", "tf", "linsight", "gpn_msa", "cadd"]
TISS = [c for c in CACHE_COLS if c.startswith("re2g_") and c != "re2g_max"]; require(len(TISS) == 10, "tissue cols")
_num = re.compile(r"^-?[0-9.]+([eE][-+]?[0-9]+)?$")
def parse_gh(txt, gene_syms):
    if not txt: return ("", "", "")
    m = re.search(r"Name=([0-9.]+)", txt); elem = m.group(1) if m else ""
    links = re.findall(r"connected_gene=([^;]+);score=([0-9.]+)", txt)
    sc = [float(s) for g, s in links if g in gene_syms]
    return (elem, str(len(links)), (f"{max(sc)}" if sc else ""))
z = np.load(f"{A}/cache/fm_all.npz"); ours = set(z["key37"].tolist()); log(f"our cache distinct key37 {len(ours):,}"); del z
sym = {}
with gzip.open(f"{R}/deductive/gencode.sorted.gtf.gz", "rt") as fh:
    for line in fh:
        if line[0] == "#": continue
        f = line.split("\t")
        if f[2] != "gene": continue
        i = f[8].find('gene_id "'); gid = f[8][i+9:f[8].find('"', i+9)].split(".")[0]
        j = f[8].find('gene_name "'); sym[gid] = f[8][j+11:f[8].find('"', j+11)] if j >= 0 else ""
log(f"gene symbols {len(sym):,}")
st = dict(chrs=CHRS, rows=0, dropped_chr2122=0, dropped_bad_maf=0, chr38_ne_chr37=0, in_our_band=0,
          t1_na=0, t1_partial_blanked=0, no_favor=0, no_gpn=0, no_cadd=0, band=dict(ge5=0, b1_5=0, b01_1=0, lt01=0), na=dict((c, 0) for c in COLS), per_bucket={})
OUT = f"{OUTDIR}/{args.out}"; os.makedirs(OUTDIR, exist_ok=True)
with gzip.open(OUT + ".tmp", "wt") as o:
    o.write("\t".join(["variant_hg19", "chr", "pos_hg19", "pos_hg38"] + COLS + ["in_our_band", "maf_bbj", "rsq_bbj"]) + "\n")
    for M in CHRS:
        for p in (f"{A}/map38_trA/chr{M}.map.tsv", f"{A}/extract_trA/chr{M}.annot.tsv", f"{A}/extract_trA/chr{M}.cadd.tsv", f"{A}/gpn_trA/chr{M}.gpn.tsv", f"{A}/t2_trA/chr{M}.t2.tsv"): done(p)
        def load(path, exp_hdr):
            d = {}
            with open(path) as f:
                hdr = next(f).rstrip("\n").split("\t"); require(hdr == exp_hdr, f"{path} header {hdr}")
                for line in f:
                    a = line.rstrip("\n").split("\t"); require(len(a) == len(hdr), f"{path} width"); require(a[0] not in d, f"dup key {path}"); d[a[0]] = a[1:]
            return d
        fav = load(f"{A}/extract_trA/chr{M}.annot.tsv", ["key37", "cons", "epi_active", "epi_repr", "epi_trans", "tf", "cage_prom", "genehancer", "linsight"])
        cadd = {}
        with open(f"{A}/extract_trA/chr{M}.cadd.tsv") as f:
            for line in f:
                a = line.rstrip("\n").split("\t"); require(len(a) == 2 and a[0] not in cadd, "cadd row"); cadd[a[0]] = a[1]
        gpn = load(f"{A}/gpn_trA/chr{M}.gpn.tsv", ["key37", "gpn_msa"])
        t2h = ["key37", "in_body", "in_tss3kb", "in_re2g", "dist_tss", "re2g_max", "n_genes"] + TISS + ["genes"]
        t2 = load(f"{A}/t2_trA/chr{M}.t2.tsv", t2h); t2i = {c: i for i, c in enumerate(t2h[1:])}
        n_b = 0
        with open(f"{A}/map38_trA/chr{M}.t4.tsv") as f:
            hdr = next(f).rstrip("\n").split("\t"); require(hdr == ["key37_bbj", "id37", "chr37", "pos37", "chr38", "pos38", "maf", "r2", "avg_cs", "is_typed", "is_indel", "flipped"], "t4 header")
            for line in f:
                kb, id37, c37, p37, c38, p38, maf, r2, acs, ityp, iind, fl = line.rstrip("\n").split("\t")
                if c37 in ("21", "22"): st["dropped_chr2122"] += 1; continue
                try: mb = float(maf)
                except ValueError: st["dropped_bad_maf"] += 1; continue
                if not (0 <= mb <= 0.5): st["dropped_bad_maf"] += 1; continue
                if c38 != c37: st["chr38_ne_chr37"] += 1
                fv = fav.get(kb); gv = gpn.get(kb); cv = cadd.get(kb, "")
                if fv is None: st["no_favor"] += 1
                if gv is None: st["no_gpn"] += 1
                if cv == "": st["no_cadd"] += 1
                fvd = dict(zip(["cons", "epi_active", "epi_repr", "epi_trans", "tf", "cage_prom", "genehancer", "linsight"], fv)) if fv else {}
                tv = t2.get(kb); require(tv is not None, f"t2 missing for a variant in chr{M}")
                genes = tv[t2i["genes"]].split(";") if tv[t2i["genes"]] else []
                gsyms = {sym.get(g, "") for g in genes} - {""}
                gh = parse_gh(fvd.get("genehancer", ""), gsyms)
                t1 = [fvd.get("cons", ""), fvd.get("epi_active", ""), fvd.get("epi_repr", ""), fvd.get("epi_trans", ""), fvd.get("tf", ""), fvd.get("linsight", ""), (gv[0] if gv else ""), cv]
                nmiss = sum(1 for x in t1 if x == "")
                t1na = 1 if nmiss > 0 else 0
                if t1na:
                    st["t1_na"] += 1
                    if nmiss < 8: st["t1_partial_blanked"] += 1
                    t1 = [""] * 8
                for x in t1: require(x == "" or _num.match(x), "non-numeric T1")
                row = dict(zip(T1, t1))
                row.update(is_cage_prom="1" if fvd.get("cage_prom") else "0", gh_elem_score=gh[0], gh_n_genes=gh[1], gh_link_score=gh[2],
                           in_body=tv[t2i["in_body"]], in_tss3kb=tv[t2i["in_tss3kb"]], in_re2g=tv[t2i["in_re2g"]], dist_tss=tv[t2i["dist_tss"]],
                           re2g_max=tv[t2i["re2g_max"]], n_genes=tv[t2i["n_genes"]], t1_na=str(t1na), maf=maf, r2=r2, avg_cs=acs, is_typed=ityp, is_indel=iind)
                for t in TISS: row[t] = tv[t2i[t]]
                vals = [row[c] for c in COLS]
                for c, x in zip(COLS, vals):
                    if x == "": st["na"][c] += 1
                ib = 1 if id37 in ours else 0; st["in_our_band"] += ib
                if mb >= .05: st["band"]["ge5"] += 1
                elif mb >= .01: st["band"]["b1_5"] += 1
                elif mb >= .001: st["band"]["b01_1"] += 1
                else: st["band"]["lt01"] += 1
                o.write("\t".join([kb, c37, p37, p38] + vals + [str(ib), "", ""]) + "\n")
                st["rows"] += 1; n_b += 1
        st["per_bucket"][M] = n_b
        log(f"chr{M}: {n_b:,} rows (cum {st['rows']:,})")
require(st["rows"] > 0, "0 rows")
st["na_pct"] = {c: round(st["na"][c] / st["rows"] * 100, 2) for c in COLS}
os.replace(OUT + ".tmp", OUT)
json.dump(st, open(OUT + ".stats.json", "w"), indent=1)
open(OUT + ".done", "w").write(f"ok {st['rows']}\n")
log(json.dumps({k: v for k, v in st.items() if k not in ("na", "per_bucket", "na_pct")}))
print("ASSEMBLE_DONE")
