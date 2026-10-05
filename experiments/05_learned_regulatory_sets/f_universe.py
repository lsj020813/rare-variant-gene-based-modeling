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


import gzip, json, os, subprocess, sys, collections
CHR = sys.argv[1]
ROOT = _config_path("${PROJECT_ROOT}/work")
OUT = f"{ROOT}/fset/out"
BCF = "bcftools"
VP = f"{ROOT}/gate1/out/vpos_chr{CHR}.tsv.gz"
SRC = f"{ROOT}/ref/orig_index/chr{CHR}.vcf.gz"
os.makedirs(OUT, exist_ok=True)

pos38 = {}
coding = {}
with gzip.open(VP, "rt") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) != 3 or p[0] == "hg19_id":
            continue
        pos38[p[0]] = int(p[1]); coding[p[0]] = int(p[2])

def maf_bin(m):
    return "B1" if m < 0.001 else "B2" if m < 0.01 else "B3" if m < 0.05 else "B4"
def r2_band(r):
    return "lt0.3" if r < 0.3 else "0.3-0.6" if r < 0.6 else "0.6-0.8" if r < 0.8 else "0.8-0.9" if r < 0.9 else "ge0.9"

agg = collections.Counter()
n_src = n_hit = 0
tmp = f"{OUT}/uni_chr{CHR}.tsv.gz.tmp"
proc = subprocess.Popen([BCF, "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t%INFO/MAF\t%INFO/R2\t%INFO/TYPED\n", SRC],
                        stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
with gzip.open(tmp, "wt") as out:
    out.write("key\tpos38\tmaf\tr2\ttyped\tis_coding\tvar_dosage_2pq\n")
    for line in proc.stdout:
        ch, po, ref, alt, mafs, r2s, ty = line.rstrip("\n").split("\t")
        n_src += 1
        key = f"{ch}:{po}:{ref}:{alt}"
        p38 = pos38.get(key)
        if p38 is None:
            continue
        n_hit += 1
        try: m = float(mafs)
        except ValueError: m = 0.0
        try: r = float(r2s)
        except ValueError: r = 0.0
        t = 1 if ty == "1" else 0
        cd = coding.get(key, 0)
        agg[(maf_bin(m), r2_band(r), "TYPED" if t else "IMPUTED", "coding" if cd else "noncoding")] += 1
        out.write(f"{key}\t{p38}\t{m:.6g}\t{r:.4g}\t{t}\t{cd}\t{2*m*(1-m)*r:.6g}\n")
proc.wait()
os.replace(tmp, f"{OUT}/uni_chr{CHR}.tsv.gz")
with open(f"{OUT}/uni_chr{CHR}.counts.tsv", "w") as f:
    f.write("chr\tmaf_bin\tr2_band\ttyped\tcoding\tn\n")
    for (b, rb, t, cd), n in sorted(agg.items()):
        f.write(f"{CHR}\t{b}\t{rb}\t{t}\t{cd}\t{n}\n")
print(json.dumps({"chr": CHR, "src_records": n_src, "assigned_matched": n_hit,
                  "assigned_expected": len(pos38), "noncoding": sum(v for k, v in agg.items() if k[3] == "noncoding")}))
