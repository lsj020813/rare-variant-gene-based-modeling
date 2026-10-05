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


import sys, os, gzip, json, bisect
from collections import defaultdict
CH = sys.argv[1]
R = _config_path("${PROJECT_ROOT}/work/ref")
GTF = f"{R}/deductive/gencode.sorted.gtf.gz"
T4 = f"{R}/annot/map38_bbj/chr{CH}.t4.tsv"
RE2G = f"{R}/re2g"
OUT = f"{R}/annot/t2_bbj"; os.makedirs(OUT, exist_ok=True); O = f"{OUT}/chr{CH}.t2.tsv"
if os.path.exists(O + ".done"): print(f"[chr{CH}] 이미 완료"); sys.exit(0)
assert os.path.exists(f"{R}/annot/map38_bbj/chr{CH}.map.tsv.done"), "GATE FAIL: map not done"

genes = {}
with gzip.open(GTF, "rt") as fh:
    for line in fh:
        if line[0] == "#": continue
        f = line.split("\t")
        if f[2] != "gene" or f[0] != f"chr{CH}" or 'gene_type "protein_coding"' not in f[8]: continue
        i = f[8].find('gene_id "'); gid = f[8][i+9:f[8].find('"', i+9)]
        s, e = int(f[3]), int(f[4])
        genes[gid] = {"tss": s if f[6] == "+" else e, "start": s, "end": e}
assert genes, "GATE FAIL: no protein_coding genes"
base2full = {g.split(".")[0]: g for g in genes}

var = []
with open(T4) as fh:
    hdr = next(fh).rstrip("\n").split("\t"); ik, ip, ic = hdr.index("key37_bbj"), hdr.index("pos38"), hdr.index("chr38")
    for line in fh:
        a = line.rstrip("\n").split("\t"); assert a[ic] == CH, "GATE FAIL: chr38 mismatch in t4 bucket"
        var.append((a[ik], int(a[ip])))
assert var, "GATE FAIL: 0 variants"

links = defaultdict(list); tissues = []
BIN = 10000; binidx = defaultdict(list)
for fp in sorted(os.listdir(RE2G)):
    if not fp.endswith(".bed.gz"): continue
    with gzip.open(f"{RE2G}/{fp}", "rt") as fh:
        hdr = next(fh).rstrip("\n").split("\t")
        assert hdr[13] == "TargetGeneEnsembl_ID" and hdr[27] == "CellType" and hdr[29] == "Score", f"GATE FAIL: header {fp}"
        ct = None
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 30 or f[0] != f"chr{CH}": continue
            g = base2full.get(f[13].split(".")[0])
            if not g: continue
            ct = f[27]; es, ee, sc = int(f[1]), int(f[2]), float(f[29])
            links[g].append((es, ee, ct, sc))
            for b in range(es // BIN, ee // BIN + 1): binidx[b].append((es, ee, g, ct, sc))
        if ct: tissues.append(ct)
assert len(tissues) == 10 and len(set(tissues)) == 10, f"GATE FAIL: tissues {tissues}"

gl = sorted(genes.items(), key=lambda kv: kv[1]["start"])
starts = [gi["start"] for _, gi in gl]
maxlen = max(gi["end"] - gi["start"] for _, gi in gl)
n_assigned = 0; n_zero = 0; dist_src = defaultdict(int)
with open(O + ".tmp", "w") as o:
    o.write("key37\tin_body\tin_tss3kb\tin_re2g\tdist_tss\tre2g_max\tn_genes\t" + "\t".join("re2g_" + t.replace(" ", "_") for t in tissues) + "\tgenes\n")
    for v, p in var:
        assigned = set(); body = tss = 0; dists = []
        hi = bisect.bisect_right(starts, p + 3000); lo = bisect.bisect_left(starts, p - 3000 - maxlen)
        for g, gi in gl[lo:hi]:
            b = gi["start"] <= p <= gi["end"]; t = abs(p - gi["tss"]) <= 3000
            if b or t:
                assigned.add(g); body |= b; tss |= t
        per = {}; regenes = set()
        for es, ee, g, ct, sc in binidx.get(p // BIN, ()):
            if es <= p <= ee:
                per[ct] = max(per.get(ct, -1.0), sc); regenes.add(g)
        assigned |= regenes
        re = int(bool(per))
        for g in assigned: dists.append(abs(p - genes[g]["tss"]))
        if assigned: n_assigned += 1
        else: n_zero += 1
        mx = f"{max(per.values()):.4f}" if per else ""
        o.write(f"{v}\t{int(body)}\t{int(tss)}\t{re}\t{min(dists) if dists else ''}\t{mx}\t{len(assigned)}\t" +
                "\t".join(f"{per[t]:.4f}" if t in per else "" for t in tissues) + "\t" + ";".join(sorted(g.split('.')[0] for g in assigned)) + "\n")
os.replace(O + ".tmp", O)
open(f"{O}.done", "w").write(f"ok {len(var)}\n")
print("T2STATS", json.dumps({"chr": CH, "variants": len(var), "assigned": n_assigned, "unassigned": n_zero, "tissues": tissues, "genes": len(genes)}))
print(f"[chr{CH}] T2_DONE")
