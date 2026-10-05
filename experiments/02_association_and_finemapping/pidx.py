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


import gzip, os, sys, json, glob, subprocess
from multiprocessing import Pool

THR = 0.1
OUT   = _config_path("${PROJECT_ROOT}/work/ref/gate4/pangolin_idx")
KEYED = _config_path("${PROJECT_ROOT}/work/ref/lift38_keyed")
PANG  = _config_path("${PROJECT_ROOT}/work/ref/pangolin/Pangolin_hg38_snvs_masked")
BCF   = "bcftools"
os.makedirs(OUT, exist_ok=True)
SMOKE = "--smoke" in sys.argv

def band_keys():
    ks = set()
    chroms = ["21"] if SMOKE else [str(i) for i in range(1, 23)]
    for n in chroms:
        p = subprocess.Popen(
            f"{BCF} query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' {KEYED}/chr{n}.keyed38.vcf.gz",
            shell=True, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
        for line in p.stdout:
            f = line.rstrip("\n").split("\t")
            if len(f) == 4:
                ks.add((f[0].replace("chr", ""), f[1], f[2], f[3]))
        p.wait()
    return ks

KEYS = None
def init(k):
    global KEYS
    KEYS = k

def one(fp):
    hits = []
    try:
        with gzip.open(fp, "rt") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            ix = {c: i for i, c in enumerate(hdr)}
            if any(c not in ix for c in ("chrom","pos","ref","alt","gain_score","loss_score")):
                return hits
            for line in fh:
                f = line.rstrip("\n").split("\t")
                try:
                    ch = f[ix["chrom"]].replace("chr", "")
                    k = (ch, f[ix["pos"]], f[ix["ref"]], f[ix["alt"]])
                    if k not in KEYS: continue
                    a, b = abs(float(f[ix["gain_score"]])), abs(float(f[ix["loss_score"]]))
                except (ValueError, IndexError):
                    continue
                v, d = (a, 1) if a >= b else (b, -1)
                if v >= THR:
                    hits.append((ch, ":".join(k), d, v))
    except (OSError, EOFError):
        pass
    return hits

if __name__ == "__main__":
    keys = band_keys()
    print(f"band keys: {len(keys):,}", flush=True)
    files = sorted(glob.glob(f"{PANG}/*.tsv.gz"))
    if SMOKE:
        want, picked = 25, []
        for fp in files:
            try:
                with gzip.open(fp, "rt") as fh:
                    fh.readline()
                    row = fh.readline().split("\t")
                if row and row[0].replace("chr", "") == "21":
                    picked.append(fp)
                    if len(picked) >= want: break
            except (OSError, EOFError):
                continue
        files = picked
        print(f"smoke: {len(files)} chr21 gene files selected", flush=True)
    print(f"pangolin files: {len(files):,}", flush=True)
    per = {}
    with Pool(6, initializer=init, initargs=(keys,)) as pool:
        for i, hits in enumerate(pool.imap_unordered(one, files, chunksize=20), 1):
            for ch, k, d, v in hits:
                per.setdefault(ch, []).append(f"{k}\t{d}\t{v}")
            if i % 2000 == 0: print(f"  {i:,}/{len(files):,}", flush=True)
    tot = 0
    suffix = ".smoke" if SMOKE else ""
    for ch, rows in per.items():
        with open(f"{OUT}/chr{ch}{suffix}.tsv", "w") as fh:
            fh.write("\n".join(sorted(set(rows))) + "\n")
        tot += len(set(rows))
    print(f"indexed variants: {tot:,} across {len(per)} chroms")
    if tot == 0:
        print("GATE FAIL: zero Pangolin/band intersections — key format mismatch"); sys.exit(6)
    json.dump({"variants": tot, "chroms": sorted(per), "files": len(files), "smoke": SMOKE},
              open(f"{OUT}/summary{suffix}.json", "w"), indent=1)
    print("PANGOLIN_INDEX_OK")
