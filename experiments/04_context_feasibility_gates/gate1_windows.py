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


import gzip, json, os, re, sys, subprocess, collections

CHR = sys.argv[1]
ROOT = _config_path("${PROJECT_ROOT}/work")
OUT = f"{ROOT}/gate1/out"
GENCODE = f"{ROOT}/ref/deductive/gencode.nochr.gtf.gz"
RE2G_DIR = f"{ROOT}/ref/re2g"
GROUPF = f"{ROOT}/ref/groupfiles_bwg/chr{CHR}.B_3kb_re2g.txt"
KEYED = f"{ROOT}/ref/lift38_keyed/chr{CHR}.keyed38.vcf.gz"
BCF = "bcftools"
FLANK = 3000

def bin_of(maf):
    m = float(maf)
    if m < 0.001: return "B1"
    if m < 0.01:  return "B2"
    if m < 0.05:  return "B3"
    return "B4"

genes = []
with open(GROUPF) as f:
    for line in f:
        p = line.split()
        if len(p) > 2 and p[1] == "var":
            genes.append(p[0])
genes = sorted(set(genes))
gset_nov = {g.split(".")[0]: g for g in genes}

body = {}
with gzip.open(GENCODE, "rt") as f:
    for line in f:
        if line[0] == "#": continue
        p = line.rstrip("\n").split("\t")
        if p[2] != "gene" or p[0] != CHR: continue
        m = re.search(r'gene_id "([^"]+)"', p[8])
        if not m: continue
        gid = m.group(1); nov = gid.split(".")[0]
        if nov in gset_nov:
            s, e = int(p[3]), int(p[4])
            key = gset_nov[nov]
            if key in body:
                s = min(s, body[key][0]); e = max(e, body[key][1])
            body[key] = (s, e)

re2g = collections.defaultdict(list)
for fn in sorted(os.listdir(RE2G_DIR)):
    if not fn.endswith(".bed.gz"): continue
    with gzip.open(os.path.join(RE2G_DIR, fn), "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        try:
            i_c, i_s, i_e = 0, 1, 2
            i_g = hdr.index("TargetGeneEnsembl_ID"); i_t = hdr.index("CellType")
        except ValueError:
            continue
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[i_c] not in (f"chr{CHR}", CHR): continue
            nov = p[i_g].split(".")[0]
            if nov not in gset_nov: continue
            re2g[gset_nov[nov]].append((int(p[i_s]), int(p[i_e]), p[i_t]))

def merge(iv):
    iv = sorted(iv)
    out = []
    for s, e, src in iv:
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
            out[-1][2] = out[-1][2] | src
        else:
            out.append([s, e, src])
    return out

win = {}
for g in genes:
    iv = []
    if g in body:
        iv.append((body[g][0] - FLANK, body[g][1] + FLANK, 1))
    for s, e, t in re2g.get(g, []):
        iv.append((s, e, 2))
    if iv:
        win[g] = merge(iv)

flat = []
for g, ivs in win.items():
    for s, e, src in ivs:
        flat.append((s, e, g, src))
flat.sort()
starts = [x[0] for x in flat]
import bisect

assign = collections.defaultdict(lambda: collections.defaultdict(list))
n_seen = n_assigned = 0
proc = subprocess.Popen([BCF, "query", "-f", "%CHROM\t%POS\t%ID\t%INFO/MAF\t%INFO/R2\n", KEYED],
                        stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
for line in proc.stdout:
    c, pos, vid, maf, r2 = line.rstrip("\n").split("\t")
    if c not in (f"chr{CHR}", CHR):
        continue
    n_seen += 1
    p = int(pos)
    i = bisect.bisect_right(starts, p) - 1
    j = i; hits = {}
    while j >= 0 and j > i - 400:
        s, e, g, src = flat[j]
        if e >= p and s <= p:
            hits[g] = hits.get(g, 0) | src
        j -= 1
    if not hits: continue
    n_assigned += 1
    b = bin_of(maf)
    for g, src in hits.items():
        assign[g][b].append((vid, float(maf), float(r2), src))
proc.wait()

orig_n = int(subprocess.run([BCF, "index", "-n", f"{ROOT}/ref/orig_index/chr{CHR}.vcf.gz"],
                            capture_output=True, text=True).stdout.strip() or 0)
cov = n_seen / orig_n if orig_n else 0.0
if cov < 0.95:
    raise SystemExit(f"FAIL keyed38 coverage chr{CHR}: on-chr sites {n_seen} vs orig_index {orig_n} = {cov:.3f} < 0.95")

os.makedirs(OUT, exist_ok=True)
rows = []
payload = {}
for g in genes:
    d = assign.get(g, {})
    ent = {"span_hg38": win.get(g, []), "bins": {}}
    for b in ("B1", "B2", "B3", "B4"):
        lst = d.get(b, [])
        ent["bins"][b] = [[v, m, r, s] for v, m, r, s in lst]
    payload[g] = ent
    rows.append([g, len(d.get("B1", [])), len(d.get("B2", [])), len(d.get("B3", [])), len(d.get("B4", [])),
                 sum(1 for _, _, _, s in d.get("B2", []) if s == 1),
                 sum(1 for _, _, _, s in d.get("B2", []) if s == 2),
                 sum(1 for _, _, _, s in d.get("B2", []) if s == 3),
                 sum(e - s for s, e, _ in win.get(g, [])), len(win.get(g, []))])

with gzip.open(f"{OUT}/windows_chr{CHR}.json.gz", "wt") as f:
    json.dump(payload, f)
with open(f"{OUT}/window_summary_chr{CHR}.tsv.tmp", "w") as f:
    f.write("gene\tnB1\tnB2\tnB3\tnB4\tB2_dist_only\tB2_re2g_only\tB2_both\twindow_bp\tn_intervals\n")
    for r in rows:
        f.write("\t".join(map(str, r)) + "\n")
os.replace(f"{OUT}/window_summary_chr{CHR}.tsv.tmp", f"{OUT}/window_summary_chr{CHR}.tsv")
print(json.dumps({"chr": CHR, "genes": len(genes), "genes_with_window": len(win),
                  "sites_on_chr": n_seen, "orig_index_n": orig_n, "keyed38_coverage": round(cov, 4),
                  "sites_assigned": n_assigned,
                  "genes_with_re2g": sum(1 for g in genes if re2g.get(g))}))
