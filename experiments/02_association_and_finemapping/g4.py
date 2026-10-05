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


import sys, gzip, os, json, subprocess, glob

SMOKE = "--smoke" in sys.argv
CH = ["21"] if SMOKE else [str(i) for i in range(1, 23)]
THR = 0.1
OUT = _config_path("${PROJECT_ROOT}/work/ref/gate4")
BAND = _config_path("${PROJECT_ROOT}/work/ref/band_vcf")
KEYED = _config_path("${PROJECT_ROOT}/work/ref/lift38_keyed")
SAI = (_config_path("${PROJECT_ROOT}/"
       "resources/annotation_data/spliceAI/spliceai_scores.raw.snv.hg38.vcf.gz"))
PANG_DIR = _config_path("${PROJECT_ROOT}/work/ref/pangolin/Pangolin_hg38_snvs_masked")
BCF = "bcftools"
os.makedirs(OUT, exist_ok=True)

def sh(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout

def band_keys(n):
    txt = sh(f"{BCF} query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' {KEYED}/chr{n}.keyed38.vcf.gz")
    s = set()
    for line in txt.splitlines():
        f = line.split("\t")
        if len(f) == 4:
            s.add((f[0].replace("chr", ""), int(f[1]), f[2], f[3]))
    return s

def spliceai_dir(n, keys):
    out = {}
    cmd = (f"{BCF} query -r {n} -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/SpliceAI\\n' {SAI} 2>/dev/null")
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True, bufsize=1<<20)
    for line in p.stdout:
        f = line.rstrip("\n").split("\t")
        if len(f) < 5: continue
        k = (f[0].replace("chr",""), int(f[1]), f[2], f[3])
        if k not in keys: continue
        best_v, best_d = 0.0, 0
        for ann in f[4].split(","):
            p2 = ann.split("|")
            if len(p2) < 6: continue
            try: ag, al, dg, dl = (float(x) if x not in ("", ".") else 0.0 for x in p2[2:6])
            except ValueError: continue
            for v, d in ((ag, +1), (dg, +1), (al, -1), (dl, -1)):
                if v > best_v: best_v, best_d = v, d
        if best_v >= THR:
            out[k] = (best_d, best_v)
    p.wait()
    return out

def pangolin_dir(n, keys):
    out = {}
    for fp in glob.glob(f"{PANG_DIR}/*.tsv.gz"):
        try:
            with gzip.open(fp, "rt") as fh:
                hdr = fh.readline().rstrip("\n").split("\t")
                idx = {c: i for i, c in enumerate(hdr)}
                need = ("chrom","pos","ref","alt","gain_score","loss_score")
                if any(c not in idx for c in need): continue
                inc, dec = "gain_score", "loss_score"
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if f[idx["chrom"]].replace("chr","") != n: continue
                    try: k = (n, int(f[idx["pos"]]), f[idx["ref"]], f[idx["alt"]])
                    except (ValueError, IndexError): continue
                    if k not in keys: continue
                    try: a, b = abs(float(f[idx[inc]])), abs(float(f[idx[dec]]))
                    except (ValueError, IndexError): continue
                    v, d = (a, +1) if a >= b else (b, -1)
                    if v >= THR: out[k] = (d, v)
        except (OSError, EOFError):
            continue
    return out

res = {}
for n in CH:
    keys = band_keys(n)
    sa = spliceai_dir(n, keys)
    pg = pangolin_dir(n, keys)
    both = set(sa) & set(pg)
    agree = sum(1 for k in both if sa[k][0] == pg[k][0])
    res[f"chr{n}"] = {"band": len(keys), "spliceai_dir": len(sa), "pangolin_dir": len(pg),
                      "both": len(both), "agree": agree}
    print(f"[chr{n}] band {len(keys):,} | SpliceAI {len(sa):,} | Pangolin {len(pg):,} "
          f"| both {len(both):,} | agree {agree:,}", flush=True)

T = {k: sum(v[k] for v in res.values()) for k in ("band","spliceai_dir","pangolin_dir","both","agree")}
res["TOTAL"] = T
json.dump(res, open(f"{OUT}/gate4{'_smoke' if SMOKE else ''}.json","w"), indent=1)

fail = []
if T["spliceai_dir"] == 0: fail.append("SpliceAI produced ZERO directional calls — key/format mismatch")
if T["pangolin_dir"] == 0: fail.append("Pangolin produced ZERO directional calls — key/format mismatch")
if fail:
    for m in fail: print("GATE FAIL:", m)
    sys.exit(6)
rate = T["agree"]/T["both"] if T["both"] else float("nan")
print(f"\nTOTAL band {T['band']:,} | both-scored {T['both']:,} | agree {T['agree']:,} | rate {rate:.4f}")
print("SMOKE_OK" if SMOKE else "GATE4_COMPLETE")
