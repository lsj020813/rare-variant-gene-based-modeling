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


import sys, os, gzip, json, bisect, subprocess
from collections import defaultdict
CH = sys.argv[1]
R = _config_path("${PROJECT_ROOT}/work/ref15")
GTF = f"{R}/deductive/gencode.sorted.gtf.gz"
KEYED = f"{R}/lift38_keyed/chr{CH}.keyed38.vcf.gz"
GP = f"{R}/groupfiles_bwg/chr{CH}.B_3kb_re2g.txt"
RE2G = f"{R}/re2g"
BCF = "bcftools"
OUT = f"{R}/annot/t2"; os.makedirs(OUT, exist_ok=True); O = f"{OUT}/chr{CH}.t2.tsv"

genes = {}
with gzip.open(GTF, "rt") as fh:
    for line in fh:
        if line[0] == "#": continue
        f = line.split("\t")
        if f[2] != "gene" or f[0] != f"chr{CH}" or 'gene_type "protein_coding"' not in f[8]: continue
        i = f[8].find('gene_id "'); gid = f[8][i+9:f[8].find('"', i+9)]
        s, e = int(f[3]), int(f[4])
        genes[gid] = {"tss": s if f[6] == "+" else e, "start": s, "end": e}
base2full = {g.split(".")[0]: g for g in genes}

pos38 = {}
txt = subprocess.run(f"{BCF} query -f '%ID\t%CHROM\t%POS\n' {KEYED}", shell=True, capture_output=True, text=True).stdout
for line in txt.splitlines():
    f = line.split("\t")
    if len(f) != 3 or f[1].replace("chr", "") != CH: continue
    k = f[0][3:] if f[0].startswith("chr") else f[0]
    pos38[k] = int(f[2])

links = defaultdict(list); tissues = []
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
            ct = f[27]
            links[g].append((int(f[1]), int(f[2]), ct, float(f[29])))
        if ct: tissues.append(ct)
assert len(tissues) == 10 and len(set(tissues)) == 10, f"GATE FAIL: tissues {tissues}"

pairs = []
with open(GP) as f:
    for line in f:
        a = line.split()
        if len(a) >= 3 and a[1] == "var":
            for v in a[2:]: pairs.append((v, a[0]))
ngenes = defaultdict(int)
for v, g in pairs: ngenes[v] += 1

unassigned = 0
with open(O, "w") as o:
    o.write("key37\tgene\tin_body\tin_tss3kb\tin_re2g\tdist_tss\tre2g_max\tn_genes\t" + "\t".join("re2g_" + t.replace(" ", "_") for t in tissues) + "\n")
    for v, g in pairs:
        gi = genes[g]; p = pos38[v]
        body = int(gi["start"] <= p <= gi["end"]); tss = int(abs(p - gi["tss"]) <= 3000)
        per = {}
        for es, ee, ct, sc in links.get(g, ()):
            if es <= p <= ee: per[ct] = max(per.get(ct, -1.0), sc)
        re = int(bool(per))
        if not (body or tss or re): unassigned += 1
        mx = f"{max(per.values()):.4f}" if per else ""
        o.write(f"{v}\t{g}\t{body}\t{tss}\t{re}\t{abs(p-gi['tss'])}\t{mx}\t{ngenes[v]}\t" +
                "\t".join(f"{per[t]:.4f}" if t in per else "" for t in tissues) + "\n")
assert unassigned == 0, f"GATE FAIL: {unassigned} pairs match no source (bwg rule not reproduced)"
open(f"{O}.done", "w").write(f"ok {len(pairs)}\n")
print("T2STATS", json.dumps({"chr": CH, "pairs": len(pairs), "tissues": tissues}))
print(f"[chr{CH}] T2_DONE")
