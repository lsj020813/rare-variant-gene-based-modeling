#!/usr/bin/env python3
import os as _cfg_os
import re as _cfg_re
import math as _cfg_math

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)

def _config_number(name, cast=int, positive=False):
    value = cast(_cfg_os.environ[name])
    if not _cfg_math.isfinite(value) or (positive and value <= 0):
        raise ValueError(name + " must be finite" + (" and positive" if positive else ""))
    return value

import glob, gzip, os, zipfile, json, collections, sys, re
D = _config_path("${PROJECT_ROOT}/work/ref/bbj_fm"); OUT = _config_path("${PROJECT_ROOT}/work/ref/annot/bbj_l1/raw"); os.makedirs(f"{OUT}/pip_raw", exist_ok=True)
zips = sorted(glob.glob(f"{D}/hum0197.v5.finemap.*.v1.zip"))
af = {}; afn = collections.Counter(); afdiff = collections.defaultdict(float); afsum = collections.defaultdict(float)
dup = {}; tr_rows = []
NEED = ["chromosome","position","allele1","allele2","af_allele2","beta_marginal","se_marginal","region","pip","cs_id"]
for z in zips:
    tr = z.split(".finemap.")[1].split(".v1")[0]
    assert re.fullmatch(r"[A-Za-z0-9_]+", tr), f"GATE FAIL: trait code {tr}"
    zf = zipfile.ZipFile(z)
    names = [n for n in zf.namelist() if "SuSiE" in n and n.endswith(".tsv.gz")]
    assert len(names) == 1, f"GATE FAIL: {tr} SuSiE files {names}"
    seen = set(); ndup = 0; n = 0; regions = set(); cs = set(); csdist = collections.Counter(); nonnum = 0; chrs = set()
    with gzip.open(zf.open(names[0]), "rt") as fh, open(f"{OUT}/pip_raw/{tr}.tsv.tmp", "w") as o:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {c: hdr.index(c) for c in NEED}
        o.write("variant_hg19\tchr\tpos\tallele1\tallele2\tpip\tcs_id\tregion\tbeta_marginal\tse\n")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            try: pip = float(f[ix["pip"]]); a = float(f[ix["af_allele2"]])
            except ValueError: nonnum += 1; continue
            ch = f[ix["chromosome"]]; ch = ch[3:] if ch.startswith("chr") else ch
            k = f"{ch}:{f[ix['position']]}:{f[ix['allele1']]}:{f[ix['allele2']]}"
            if k in seen: ndup += 1
            seen.add(k); n += 1; chrs.add(ch)
            reg = f[ix["region"]]; regions.add(reg); cid = f[ix["cs_id"]]; csdist[cid] += 1
            if cid != "-1": cs.add((reg, cid))
            o.write(f"{k}\t{ch}\t{f[ix['position']]}\t{f[ix['allele1']]}\t{f[ix['allele2']]}\t{f[ix['pip']]}\t{cid}\t{reg}\t{f[ix['beta_marginal']]}\t{f[ix['se_marginal']]}\n")
            if k in af: afdiff[k] = max(afdiff[k], abs(af[k] - a))
            else: af[k] = a
            afn[k] += 1; afsum[k] += a
    os.replace(f"{OUT}/pip_raw/{tr}.tsv.tmp", f"{OUT}/pip_raw/{tr}.tsv")
    dup[tr] = ndup
    tr_rows.append(dict(trait=tr, rows=n, variants=len(seen), regions=len(regions), cs=len(cs), nonnumeric_skipped=nonnum,
                        cs_id_values=len(csdist), cs_neg1=csdist.get("-1", 0), chrs=sorted(chrs, key=lambda x: (len(x), x))))
    print(f"[A] {tr}: 행 {n:,} 변이 {len(seen):,} 구역 {len(regions):,} CS {len(cs):,} 중복(trait,variant) {ndup} cs=-1 {csdist.get('-1',0):,} 염색체 {len(chrs)}", flush=True)
with open(f"{OUT}/variant_af.tsv.tmp", "w") as o:
    o.write("key37_bbj\taf_mean\taf_first\tn_traits\tmax_absdiff\n")
    for k, a in af.items(): o.write(f"{k}\t{afsum[k]/afn[k]:.8g}\t{a}\t{afn[k]}\t{afdiff.get(k, 0.0)}\n")
os.replace(f"{OUT}/variant_af.tsv.tmp", f"{OUT}/variant_af.tsv")
incons = sum(1 for k in af if afdiff.get(k, 0) > 1e-6)
json.dump(dict(dup_per_trait=dup, total_dup=sum(dup.values()), variants=len(af), af_inconsistent_gt1e6=incons), open(f"{OUT}/dup_report.json", "w"), indent=1)
with open(f"{OUT}/traits_raw.tsv", "w") as o:
    o.write("\t".join(tr_rows[0].keys()) + "\n")
    for r in tr_rows: o.write("\t".join(str(v) for v in r.values()) + "\n")
assert len(af) == _config_number("N_BBJ_VARIANTS", int, True), f"GATE FAIL: 고유 변이 {len(af)} != configured expected variant count"
open(f"{OUT}/A.done", "w").write("ok\n")
print(f"[A] 고유 변이 {len(af):,} · 형질 {len(zips)} · (trait,variant) 중복 합 {sum(dup.values())} · af 불일치(>1e-6) {incons}")
print("PIP_A_DONE")
