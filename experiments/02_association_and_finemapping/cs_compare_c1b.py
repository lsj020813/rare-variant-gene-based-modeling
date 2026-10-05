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


import glob, gzip, zipfile, csv, json, os, statistics as st

W = _config_path("${PROJECT_ROOT}/work")
DIRECT = {"tchl": "TC", "dm": "T2D"}
THR = 0.5

def key(ch, pos, a1, a2): return (str(ch), int(pos), frozenset((a1, a2)))

ours = {}; span = {}
for tr in DIRECT:
    ours[tr] = {}; span[tr] = {}
    for f in sorted(glob.glob(f"{W}/run_ourfm/fm/{tr}/*.pip.tsv")):
        if "mi1000" in f: continue
        lo = hi = None; ch = None; reg = None
        with open(f) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                ch = r["chr"]; reg = r["region"]; pos = int(r["pos"])
                lo = pos if lo is None else min(lo, pos); hi = pos if hi is None else max(hi, pos)
                cs = r.get("cs_id", "NA")
                if cs in ("NA", "", "-1", None): continue
                ours[tr].setdefault((reg, cs), []).append((key(ch, pos, r["a1"], r["a2"]), float(r["pip"]), float(r["maf"])))
        if reg: span[tr][reg] = (ch, lo, hi)

def load_bbj(code, ranges):
    z = f"{W}/ref/bbj_fm/hum0197.v5.finemap.{code}.v1.zip"
    zf = zipfile.ZipFile(z)
    mem = [n for n in zf.namelist() if "SuSiE" in n and n.endswith(".gz")][0]
    out = {}
    with zf.open(mem) as raw, gzip.open(raw, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ch = r["chromosome"]; rr = ranges.get(ch)
            if not rr: continue
            pos = int(r["position"])
            if not any(lo <= pos <= hi for lo, hi in rr): continue
            cs = r.get("cs_id", "-1")
            if cs in ("-1", "NA", "", None): continue
            out.setdefault((r["region"], cs), []).append(
                (key(ch, pos, r["allele1"], r["allele2"]), float(r["pip"]), float(r["af_allele2"])))
    return out

rep = {"threshold": THR, "results": {}, "mismatch_examples": []}
rows = []
for tr, code in DIRECT.items():
    ranges = {}
    for reg, (ch, lo, hi) in span[tr].items(): ranges.setdefault(ch, []).append((lo, hi))
    bbj = load_bbj(code, ranges)
    pairs = []
    for (oreg, ocs), om in ours[tr].items():
        ok = [k for k, _, _ in om]
        ol = max(om, key=lambda x: x[1])
        best = None
        for (breg, bcs), bm in bbj.items():
            bk = [k for k, _, _ in bm]
            inter = len(set(ok) & set(bk))
            if inter == 0: continue
            j = inter / len(set(ok) | set(bk))
            if best is None or j > best[1]:
                bl = max(bm, key=lambda x: x[1])
                best = ((breg, bcs), j, bl, len(bk))
        if best:
            pairs.append(dict(trait=tr, our_region=oreg, our_cs=ocs, our_n=len(ok),
                              our_maxpip=round(ol[1], 4), our_lead=ol[0], our_lead_maf=round(ol[2], 5),
                              bbj_region=best[0][0], bbj_cs=best[0][1], bbj_n=best[3],
                              bbj_maxpip=round(best[2][1], 4), bbj_lead=best[2][0], bbj_lead_af=round(best[2][2], 5),
                              jaccard=round(best[1], 4), lead_match=(ol[0] == best[2][0])))
    hi_conf = [p for p in pairs if p["our_maxpip"] >= THR and p["bbj_maxpip"] >= THR]
    lm_all = [p["lead_match"] for p in pairs]; lm_hi = [p["lead_match"] for p in hi_conf]
    rep["results"][f"{tr}~{code}"] = dict(
        pairs=len(pairs), lead_match_all=(round(sum(lm_all)/len(lm_all), 4) if lm_all else None),
        hi_conf_pairs=len(hi_conf), lead_match_hi=(round(sum(lm_hi)/len(lm_hi), 4) if lm_hi else None),
        lead_match_hi_n=sum(lm_hi),
        our_maxpip_median=round(st.median([p["our_maxpip"] for p in pairs]), 4) if pairs else None,
        bbj_maxpip_median=round(st.median([p["bbj_maxpip"] for p in pairs]), 4) if pairs else None,
        jaccard_hi_median=(round(st.median([p["jaccard"] for p in hi_conf]), 4) if hi_conf else None))
    for p in hi_conf:
        if not p["lead_match"]:
            rep["mismatch_examples"].append({k: (str(v) if isinstance(v, tuple) else v) for k, v in p.items()})
    rows += [{k: (str(v) if isinstance(v, tuple) else v) for k, v in p.items()} for p in pairs]

with open(f"{W}/run_ourfm/c1b_hiconf.tsv", "w", newline="") as fh:
    if rows:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t"); wr.writeheader(); wr.writerows(rows)
json.dump(rep, open(f"{W}/run_ourfm/c1b_hiconf.json", "w"), ensure_ascii=False, indent=1)

print(f"=== 확신 높은 부분집합 (양쪽 max PIP >= {THR})")
for k, v in rep["results"].items():
    print(f"  {k:12s} 짝 {v['pairs']:3d} (lead 일치 {v['lead_match_all']}) | "
          f"고신뢰 짝 {v['hi_conf_pairs']:3d} → **lead 일치 {v['lead_match_hi_n']}/{v['hi_conf_pairs']} = {v['lead_match_hi']}** | "
          f"자카드 중위 {v['jaccard_hi_median']} | max PIP 중위 우리 {v['our_maxpip_median']} BBJ {v['bbj_maxpip_median']}")
print(f"=== 고신뢰 불일치 예 {len(rep['mismatch_examples'])}건 (최대 6개)")
for e in rep["mismatch_examples"][:6]:
    print(f"  {e['trait']} {e['our_region']} 우리 lead {e['our_lead']} (pip {e['our_maxpip']}, maf {e['our_lead_maf']}) "
          f"vs BBJ {e['bbj_lead']} (pip {e['bbj_maxpip']}, af {e['bbj_lead_af']}) | 자카드 {e['jaccard']}")
print("C1B_DONE")
