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
PROXY  = {"htn": ["SBP", "DBP"], "lip": ["LDLC", "TG"]}

def key(ch, pos, a1, a2):
    return (str(ch), int(pos), frozenset((a1, a2)))

ours = {}
span = {}
sum_ncs = {}
for tr in list(DIRECT) + list(PROXY):
    ours[tr] = {}; span[tr] = {}; sum_ncs[tr] = 0
    for f in sorted(glob.glob(f"{W}/run_ourfm/fm/{tr}/*.summary.tsv")):
        if "mi1000" in f: continue
        with open(f) as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            for r in rd:
                try: sum_ncs[tr] += int(r["n_cs"])
                except (ValueError, KeyError): pass
    for f in sorted(glob.glob(f"{W}/run_ourfm/fm/{tr}/*.pip.tsv")):
        if "mi1000" in f: continue
        with open(f) as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            lo = hi = None; ch = None; reg = None
            for r in rd:
                ch = r["chr"]; reg = r["region"]; pos = int(r["pos"])
                lo = pos if lo is None else min(lo, pos); hi = pos if hi is None else max(hi, pos)
                cs = r.get("cs_id", "NA")
                if cs in ("NA", "", "-1", None): continue
                k = key(ch, pos, r["a1"], r["a2"])
                ours[tr].setdefault((reg, cs), []).append((k, float(r["pip"])))
            if reg: span[tr][reg] = (ch, lo, hi)

def load_bbj(trait_code, ranges):
    z = f"{W}/ref/bbj_fm/hum0197.v5.finemap.{trait_code}.v1.zip"
    if not os.path.exists(z): return None
    zf = zipfile.ZipFile(z)
    mem = [n for n in zf.namelist() if "SuSiE" in n and n.endswith(".gz")]
    if not mem: return None
    out = {}; n = 0
    with zf.open(mem[0]) as raw, gzip.open(raw, "rt") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            n += 1
            ch = r["chromosome"]
            rr = ranges.get(ch)
            if not rr: continue
            pos = int(r["position"])
            if not any(lo <= pos <= hi for lo, hi in rr): continue
            cs = r.get("cs_id", "-1")
            if cs in ("-1", "NA", "", None): continue
            out.setdefault((r["region"], cs), []).append((key(ch, pos, r["allele1"], r["allele2"]), float(r["pip"])))
    return out, n

def jaccard(a, b):
    A, B = set(a), set(b)
    return len(A & B) / len(A | B) if (A | B) else 0.0

def compare(our_cs, bbj_cs):
    rows = []; matched_our = set(); matched_bbj = set()
    for (oreg, ocs), omem in our_cs.items():
        ok = [k for k, _ in omem]
        olead = max(omem, key=lambda x: x[1])[0]
        best = None
        for (breg, bcs), bmem in bbj_cs.items():
            bk = [k for k, _ in bmem]
            inter = len(set(ok) & set(bk))
            if inter == 0: continue
            j = jaccard(ok, bk)
            if best is None or j > best[2]:
                blead = max(bmem, key=lambda x: x[1])[0]
                best = ((breg, bcs), inter, j, blead, len(bk))
        if best:
            matched_our.add((oreg, ocs)); matched_bbj.add(best[0])
            rows.append(dict(our_region=oreg, our_cs=ocs, our_n=len(ok),
                             bbj_region=best[0][0], bbj_cs=best[0][1], bbj_n=best[4],
                             inter=best[1], jaccard=round(best[2], 4),
                             lead_match=(olead == best[3])))
    return rows, matched_our, matched_bbj

report = {"trait_map": {"direct": DIRECT, "proxy": PROXY},
          "note": "BBJ SuSiE 파일 사용(우리도 SuSiE). 대립 순서 무관 매칭. pip.tsv cs_id 과소집계 한계 명시.",
          "results": {}, "our_cs_counts": {}, "summary_ncs": sum_ncs}
all_rows = []
for tr in list(DIRECT) + list(PROXY):
    ranges = {}
    for reg, (ch, lo, hi) in span[tr].items():
        ranges.setdefault(ch, []).append((lo, hi))
    report["our_cs_counts"][tr] = len(ours[tr])
    codes = [DIRECT[tr]] if tr in DIRECT else PROXY[tr]
    for code in codes:
        got = load_bbj(code, ranges)
        if got is None:
            report["results"][f"{tr}~{code}"] = {"error": "zip/member 없음"}; continue
        bbj_cs, nlines = got
        rows, mo, mb = compare(ours[tr], bbj_cs)
        for x in rows: x["pair"] = f"{tr}~{code}"
        all_rows += rows
        leads = [x["lead_match"] for x in rows]
        js = [x["jaccard"] for x in rows]
        report["results"][f"{tr}~{code}"] = dict(
            kind=("direct" if tr in DIRECT else "proxy"),
            our_cs=len(ours[tr]), bbj_cs_in_our_regions=len(bbj_cs),
            matched_pairs=len(rows), our_only=len(ours[tr]) - len(mo),
            bbj_only=len(bbj_cs) - len(mb),
            lead_match_n=sum(leads), lead_match_rate=(round(sum(leads)/len(leads), 4) if leads else None),
            jaccard_median=(round(st.median(js), 4) if js else None),
            jaccard_q1=(round(sorted(js)[len(js)//4], 4) if js else None),
            jaccard_q3=(round(sorted(js)[3*len(js)//4], 4) if js else None),
            bbj_lines_scanned=nlines)

with open(f"{W}/run_ourfm/c1_cs_compare.tsv", "w", newline="") as fh:
    if all_rows:
        wr = csv.DictWriter(fh, fieldnames=list(all_rows[0].keys()), delimiter="\t")
        wr.writeheader(); wr.writerows(all_rows)
json.dump(report, open(f"{W}/run_ourfm/c1_cs_compare.json", "w"), ensure_ascii=False, indent=1)

print("=== 우리 CS (pip.tsv cs_id 기준) vs summary n_cs 합")
for tr in report["our_cs_counts"]:
    print(f"  {tr:5s} pip.tsv {report['our_cs_counts'][tr]:4d} | summary n_cs {sum_ncs[tr]:4d}")
print("=== 비교")
for k, v in report["results"].items():
    if "error" in v: print(f"  {k}: {v['error']}"); continue
    print(f"  {k:14s} [{v['kind']:6s}] 우리 CS {v['our_cs']:3d} | BBJ CS(우리 구역내) {v['bbj_cs_in_our_regions']:4d} | 짝 {v['matched_pairs']:3d} | "
          f"lead 일치 {v['lead_match_n']:3d} ({v['lead_match_rate']}) | 자카드 중위 {v['jaccard_median']} [{v['jaccard_q1']}, {v['jaccard_q3']}] | "
          f"우리만 {v['our_only']:3d} | BBJ만 {v['bbj_only']:4d}")
print("C1_DONE")
