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


import gzip, os, sys, json, collections
ROOT = _config_path("${PROJECT_ROOT}/work"); OUT = f"{ROOT}/fset/primary"
GN = f"{ROOT}/ref/gnomad_constraint/gnomad.v4.1.constraint_metrics.tsv"
S = sys.argv[1]
by_chr = collections.defaultdict(list)
with open(S) as f:
    hdr = f.readline().rstrip("\n").split("\t"); ix = {h: i for i, h in enumerate(hdr)}
    for line in f:
        p = line.rstrip("\n").split("\t")
        by_chr[p[ix["chr"]]].append((p[ix["key"]], p[ix["pos38"]], p[ix["domain"]]))
gn = {}
with open(GN) as f:
    h = f.readline().rstrip("\n").split("\t"); gi = {c: i for i, c in enumerate(h)}
    for line in f:
        p = line.rstrip("\n").split("\t")
        g = p[gi["gene_id"]]
        if not g.startswith("ENSG"): continue
        pri = (0 if (p[gi["mane_select"]] == "true" and p[gi["canonical"]] == "true") else
               1 if p[gi["canonical"]] == "true" else 2 if p[gi["mane_select"]] == "true" else 3)
        rec = (pri, [p[gi[c]] for c in ["lof.oe", "lof.oe_ci.upper", "lof.pLI", "lof.z_score", "mis.z_score", "syn.z_score", "constraint_flags"]])
        if g not in gn or rec[0] < gn[g][0]: gn[g] = rec
def clean(v): return "" if v in ("NA", "") else v
GCOLS = ["gn_lof_oe", "gn_loeuf", "gn_pli", "gn_lof_z", "gn_mis_z", "gn_syn_z", "gn_flags"]
SCOLS = ["st_cons", "st_epi_active", "st_epi_repr", "st_epi_trans", "st_tf", "st_cage_prom", "st_genehancer", "st_linsight"]
CCOLS = ["cadd_raw", "cadd_phred"]
summary = {}
for ch in sorted(by_chr, key=lambda c: int(c) if c.isdigit() else 99):
    rows = by_chr[ch]; want = {k.replace("chr", "", 1) for k, _, _ in rows}
    st = {}
    ap = f"{ROOT}/ref/annot/extract/chr{ch}.annot.tsv"
    if os.path.exists(ap):
        with open(ap) as f:
            f.readline()
            for line in f:
                p = line.rstrip("\n").split("\t")
                if p[0] in want: st[p[0]] = [clean(x) for x in p[1:9]]
    cd = {}
    cp = f"{ROOT}/ref/features/chr{ch}.features.tsv.gz"
    if os.path.exists(cp):
        with gzip.open(cp, "rt") as f:
            h = f.readline().rstrip("\n").split("\t"); ci = {c: i for i, c in enumerate(h)}
            for line in f:
                p = line.rstrip("\n").split("\t")
                if p[0] in want: cd[p[0]] = [clean(p[ci["cadd_raw"]]), clean(p[ci["cadd_phred"]])]
    tmp = f"{OUT}/vset_chr{ch}.tsv.gz.tmp"; n_gn = n_st = n_cd = 0
    with gzip.open(tmp, "wt") as o:
        o.write("\t".join(["key", "chr", "pos38", "domain"] + GCOLS + SCOLS + CCOLS) + "\n")
        for k, pos, dom in rows:
            g = gn.get(dom.split(".")[0])
            gv = [clean(x) for x in g[1]] if g else [""] * 7
            k37 = k.replace("chr", "", 1)
            sv = st.get(k37, [""] * 8); cv = cd.get(k37, [""] * 2)
            n_gn += bool(g); n_st += k37 in st; n_cd += k37 in cd
            o.write("\t".join([k, ch, pos, dom] + gv + sv + cv) + "\n")
    os.replace(tmp, f"{OUT}/vset_chr{ch}.tsv.gz")
    summary[ch] = {"n": len(rows), "gnomad_hit": n_gn, "staar_hit": n_st, "cadd_hit": n_cd, "annot_exists": os.path.exists(ap), "cadd_exists": os.path.exists(cp)}
    print(json.dumps({ch: summary[ch]}), flush=True)
json.dump(summary, open(f"{OUT}/vset.meta.json", "w"), indent=1)
open(f"{OUT}/vset.done", "w").write("ok\n")
