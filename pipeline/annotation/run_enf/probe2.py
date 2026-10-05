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


import subprocess, sys, json, os

BCF   = "bcftools"
LIFT  = _config_path("${PROJECT_ROOT}/work/tmp/lift38")
KG    = ("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
         "1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/"
         "ALL.chr{N}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz")
OUT   = _config_path("${PROJECT_ROOT}/work/ref/enformer_probe")
os.makedirs(OUT, exist_ok=True)

def nc(c):
    return c[3:] if c.startswith("chr") else c

res = {}
for N in [int(x) for x in sys.argv[1:]] or [21, 22]:
    ours = {}
    p = subprocess.Popen(
        f"{BCF} query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/MAF\\n' {LIFT}/chr{N}.sites38.vcf.gz",
        shell=True, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for line in p.stdout:
        f = line.rstrip("\n").split("\t")
        try:
            maf = float(f[4])
        except (ValueError, IndexError):
            continue
        if 0.001 <= maf <= 0.01:
            ours[(int(f[1]), f[2], f[3])] = maf
    pos_ours = {k[0] for k in ours}
    print(f"[chr{N}] our band variants: {len(ours):,}", flush=True)

    exp0 = subprocess.run(f"{BCF} index -n '{KG.format(N=N)}'", shell=True,
                          capture_output=True, text=True, timeout=180).stdout.strip()
    exp0 = int(exp0) if exp0.isdigit() else None
    for attempt in range(1, 6):
        kg_pos, kg_key = set(), {}
        cmd = f"{BCF} query -f '%POS\\t%REF\\t%ALT\\t%INFO/AF\\n' '{KG.format(N=N)}' 2>/dev/null"
        q = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
        n_kg = 0
        for line in q.stdout:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            n_kg += 1
            pos = int(f[0]); kg_pos.add(pos)
            try:
                af = float(f[3].split(",")[0])
            except ValueError:
                af = float("nan")
            kg_key[(pos, f[1], f[2])] = af
        if exp0 is None or n_kg >= exp0:
            break
        print(f"[chr{N}] stream truncated {n_kg:,}/{exp0:,} — retry {attempt}", flush=True)
    exp = subprocess.run(f"{BCF} index -n '{KG.format(N=N)}'", shell=True,
                         capture_output=True, text=True, timeout=180).stdout.strip()
    exp = int(exp) if exp.isdigit() else None
    print(f"[chr{N}] 1000G records: {n_kg:,} (index says {exp:,})" if exp
          else f"[chr{N}] 1000G records: {n_kg:,} (index count unavailable)", flush=True)
    assert n_kg > 0, f"chr{N}: 1000G stream returned nothing"
    if exp is not None and n_kg < exp:
        raise SystemExit(f"GATE FAIL chr{N}: stream truncated {n_kg:,}/{exp:,} "
                         f"({n_kg/exp*100:.1f}%) — rerun; do NOT report this chromosome")

    a = sum(1 for k in ours if k[0] in kg_pos)
    matched = [(k, kg_key[k]) for k in ours if k in kg_key]
    b = len(matched)
    c = sum(1 for _, af in matched
            if af == af and min(af, 1 - af) >= 0.005)
    res[f"chr{N}"] = {
        "our_band": len(ours), "kg_records": n_kg,
        "pos_only": a, "pos_allele": b, "enformer_scored": c,
        "frac_pos_only": round(a / len(ours), 5),
        "frac_pos_allele": round(b / len(ours), 5),
        "frac_enformer_scored": round(c / len(ours), 5),
    }
    r = res[f"chr{N}"]
    print(f"[chr{N}] pos-only {a:,} ({r['frac_pos_only']*100:.1f}%) | "
          f"pos+allele {b:,} ({r['frac_pos_allele']*100:.1f}%) | "
          f"★enformer-scored {c:,} ({r['frac_enformer_scored']*100:.1f}%)", flush=True)

json.dump(res, open(f"{OUT}/intersection.json", "w"), indent=1)
tb = sum(v["our_band"] for v in res.values())
te = sum(v["enformer_scored"] for v in res.values())
print(f"\nTOTAL band {tb:,} | enformer-scored {te:,} ({te/tb*100:.2f}%)")
