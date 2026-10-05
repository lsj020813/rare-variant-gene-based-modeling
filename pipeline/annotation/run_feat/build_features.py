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


import os, sys, gzip, json, subprocess, collections

CHR   = sys.argv[1]
A     = _config_path("${PROJECT_ROOT}/work/ref/annotation_ref/annotation_data")
KEYED = _config_path(f'${{PROJECT_ROOT}}/work/ref/lift38_keyed/chr{CHR}.keyed38.vcf.gz')
BAND  = _config_path(f'${{PROJECT_ROOT}}/work/ref/band_vcf/chr{CHR}.band.vcf.gz')
OUT   = _config_path("${PROJECT_ROOT}/work/ref/features")
BCF   = "bcftools"
TABIX = "tabix"
os.makedirs(OUT, exist_ok=True)

def nk(key):
    return key[3:] if key.startswith("chr") else key

def nc(chrom):
    return chrom[3:] if chrom.startswith("chr") else chrom

def run(cmd):
    return subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True, bufsize=1<<20).stdout

band = set()
for line in run(f"{BCF} query -f '%CHROM:%POS:%REF:%ALT\\n' {BAND}"):
    band.add(nk(line.strip()))
print(f"[chr{CHR}] band variants: {len(band):,}", flush=True)

rows = {}
lo38 = hi38 = None
for line in run(f"{BCF} query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%ID\\n' {KEYED}"):
    c38, p38, ref, alt, k37 = line.rstrip("\n").split("\t")
    k37 = nk(k37)
    if k37 not in band:
        continue
    p = int(p38)
    rows[k37] = {"chrom38": nc(c38), "pos38": p, "ref": ref, "alt": alt}
    lo38 = p if lo38 is None or p < lo38 else lo38
    hi38 = p if hi38 is None or p > hi38 else hi38
assert rows, f"GATE FAIL chr{CHR}: 0 band variants mapped to hg38 — key convention mismatch"
print(f"[chr{CHR}] mapped to hg38: {len(rows):,} ({len(rows)/len(band)*100:.2f}%) | span {lo38:,}-{hi38:,}", flush=True)

by38 = collections.defaultdict(list)
for k, r in rows.items():
    by38[(r["chrom38"], r["pos38"], r["ref"], r["alt"])].append(k)

REG_BARE = f"{CHR}:{lo38}-{hi38}"
REG_CHR  = f"chr{CHR}:{lo38}-{hi38}"

n_cadd = 0
for f, tag in ((f"{A}/cadd/whole_genome_SNVs.tsv.gz", "snv"),
               (f"{A}/cadd/gnomad.genomes.r3.0.indel.tsv.gz", "indel")):
    if not os.path.exists(f): continue
    for line in run(f"{TABIX} {f} {REG_BARE} 2>/dev/null"):
        p = line.rstrip("\n").split("\t")
        if len(p) < 6: continue
        hit = by38.get((nc(p[0]), int(p[1]), p[2], p[3]))
        if hit:
            for k in hit:
                rows[k]["cadd_raw"], rows[k]["cadd_phred"] = p[4], p[5]
                n_cadd += 1
print(f"[chr{CHR}] CADD matched: {n_cadd:,}", flush=True)

n_spl = 0
for f in (f"{A}/spliceAI/spliceai_scores.raw.snv.hg38.vcf.gz",
          f"{A}/spliceAI/spliceai_scores.raw.indel.hg38.vcf.gz"):
    if not os.path.exists(f): continue
    for line in run(f"{BCF} query -r {REG_BARE} -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/SpliceAI\\n' {f} 2>/dev/null"):
        p = line.rstrip("\n").split("\t")
        if len(p) < 5: continue
        hit = by38.get((nc(p[0]), int(p[1]), p[2], p[3]))
        if not hit: continue
        fields = p[4].split("|")
        if len(fields) < 6: continue
        try:
            ag, al, dg, dl = (float(x) for x in fields[2:6])
        except ValueError:
            continue
        for k in hit:
            rows[k].update(spliceai_ag=ag, spliceai_al=al, spliceai_dg=dg, spliceai_dl=dl,
                           spliceai_max=max(ag, al, dg, dl),
                           spliceai_gain=max(ag, dg), spliceai_loss=max(al, dl))
            n_spl += 1
print(f"[chr{CHR}] spliceAI matched: {n_spl:,}", flush=True)

def tsv_join(path, region, cols, names, chr_prefixed=True):
    n = 0
    if not os.path.exists(path): return 0
    for line in run(f"{TABIX} {path} {region} 2>/dev/null"):
        p = line.rstrip("\n").split("\t")
        if len(p) <= max(cols): continue
        try: pos = int(p[1])
        except ValueError: continue
        hit = by38.get((nc(p[0]), pos, p[2], p[3]))
        if not hit: continue
        for k in hit:
            for c, nm in zip(cols, names):
                rows[k][nm] = p[c]
            n += 1
    return n
n_am = tsv_join(f"{A}/AlphaMissense/AlphaMissense_hg38.tsv.gz", REG_CHR, [8], ["alphamissense"])
n_pai = 0
pai = f"{A}/primateAI/PrimateAI_scores_v0.2_GRCh38_sorted.tsv.bgz"
if os.path.exists(pai):
    for reg in (REG_CHR, REG_BARE):
        n_pai = tsv_join(pai, reg, [10] if False else [4], ["primateai"])
        if n_pai: break
print(f"[chr{CHR}] AlphaMissense: {n_am:,} | primateAI: {n_pai:,}", flush=True)

FIELDS = ["cadd_raw","cadd_phred","spliceai_ag","spliceai_al","spliceai_dg","spliceai_dl",
          "spliceai_max","spliceai_gain","spliceai_loss","alphamissense","primateai"]
present = collections.Counter()
path = f"{OUT}/chr{CHR}.features.tsv.gz"
with gzip.open(path, "wt") as fh:
    fh.write("key37\tchrom38\tpos38\tref\talt\t" + "\t".join(FIELDS) + "\n")
    for k in sorted(rows, key=lambda x: rows[x]["pos38"]):
        r = rows[k]
        vals = []
        for f in FIELDS:
            v = r.get(f, "")
            if v != "": present[f] += 1
            vals.append(str(v))
        fh.write(f"{k}\t{r['chrom38']}\t{r['pos38']}\t{r['ref']}\t{r['alt']}\t" + "\t".join(vals) + "\n")

n = len(rows)
ledger = {"chrom": CHR, "band_variants": len(band), "mapped_hg38": n,
          "unmapped": len(band) - n,
          "present": dict(present),
          "missing_rate": {f: round(1 - present[f]/n, 5) for f in FIELDS}}
json.dump(ledger, open(f"{OUT}/chr{CHR}.ledger.json", "w"), indent=1)
print(f"[chr{CHR}] wrote {n:,} rows -> {path}")
for f in FIELDS:
    print(f"    {f:<16} present {present[f]:>9,}  missing {ledger['missing_rate'][f]*100:5.1f}%")
open(f"{OUT}/chr{CHR}.done", "w").close()
