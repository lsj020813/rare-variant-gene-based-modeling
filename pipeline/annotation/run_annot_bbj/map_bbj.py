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

import os, subprocess, sys, collections, json, time
R    = _config_path("${PROJECT_ROOT}/work/ref")
WANT = f"{R}/bbj_fm/union/keys37.txt"
KEY  = f"{R}/lift38_keyed"
OUT  = f"{R}/annot/map38_bbj"
BCF  = "bcftools"
SMOKE = "--smoke" in sys.argv
EXPECT_HIT = 6935025
os.makedirs(OUT, exist_ok=True)
T0 = time.time()
def log(m): print(f"[map_bbj {time.time()-T0:6.0f}s] {m}", flush=True)

want = set()
with open(WANT) as fh:
    for line in fh:
        k = line.strip()
        if k: want.add(k)
log(f"BBJ 합집합 37키: {len(want):,}")
assert len(want) == _config_number("N_BBJ_VARIANTS", int, True), f"GATE FAIL: 합집합 키 수 {len(want)} != configured expected variant count"

src = [20] if SMOKE else range(1, 23)
buckets = collections.defaultdict(list)
seen37 = set(); hit = flipped = dup37 = 0; other_contig = 0; tot_all = 0
for n in src:
    fp = f"{KEY}/chr{n}.keyed38.vcf.gz"
    assert os.path.exists(fp), f"GATE FAIL: 없음 {fp}"
    p = subprocess.Popen([BCF, "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t%ID\t%INFO/MAF\t%INFO/R2\t%INFO/AVG_CS\t%INFO/TYPED\n", fp],
                         stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    kept = tot = 0
    for line in p.stdout:
        tot += 1
        ch, pos, ref, alt, vid, maf, r2, acs, typed = line.rstrip("\n").split("\t")
        id37 = vid[3:] if vid.startswith("chr") else vid
        if id37 in want:
            kb = id37; fl = 0
        else:
            c37, p37, a, b = id37.split(":")
            kb = f"{c37}:{p37}:{b}:{a}"
            if kb not in want: continue
            fl = 1
        if kb in seen37: dup37 += 1; continue
        seen37.add(kb); hit += 1; flipped += fl; kept += 1
        cb = ch[3:] if ch.startswith("chr") else ch
        if cb not in {str(i) for i in range(1, 23)}: other_contig += 1; continue
        c37, p37 = kb.split(":")[:2]
        is_typed = 0 if typed in (".", "") else 1
        is_indel = 1 if len(ref) != len(alt) else 0
        buckets[cb].append((f"{cb}:{pos}:{ref}:{alt}", kb, id37, c37, p37, cb, pos, maf, r2, acs, str(is_typed), str(is_indel), str(fl)))
    p.stdout.close(); rc = p.wait()
    assert rc == 0, f"GATE FAIL: bcftools rc={rc} chr{n}"
    tot_all += tot
    log(f"chr{n}: 주행 {tot:,}  채택 {kept:,}  누적 히트 {hit:,}")

assert hit > 0, "GATE FAIL: 대응 0건"
if not SMOKE:
    assert hit == EXPECT_HIT, f"GATE FAIL: 히트 {hit:,} != union.log 실측 {EXPECT_HIT:,}"
dup38 = 0
for cb, rows in buckets.items():
    s = set()
    for r in rows:
        if r[0] in s: dup38 += 1
        s.add(r[0])
for cb, rows in sorted(buckets.items(), key=lambda kv: int(kv[0])):
    with open(f"{OUT}/chr{cb}.map.tsv.tmp", "w") as fm, open(f"{OUT}/chr{cb}.t4.tsv.tmp", "w") as ft:
        ft.write("key37_bbj\tid37\tchr37\tpos37\tchr38\tpos38\tmaf\tr2\tavg_cs\tis_typed\tis_indel\tflipped\n")
        for r in rows:
            fm.write(f"{r[0]}\t{r[1]}\n")
            ft.write("\t".join(r[1:]) + "\n")
    os.replace(f"{OUT}/chr{cb}.map.tsv.tmp", f"{OUT}/chr{cb}.map.tsv")
    os.replace(f"{OUT}/chr{cb}.t4.tsv.tmp", f"{OUT}/chr{cb}.t4.tsv")
    with open(f"{OUT}/chr{cb}.map.tsv.done", "w") as fh: fh.write(f"ok {len(rows)}\n")
    log(f"chr{cb}.map.tsv {len(rows):,}행")
stats = dict(smoke=SMOKE, scanned=tot_all, hit=hit, flipped=flipped, dup37_multi_record=dup37, dup38_keys=dup38,
             non_autosome38=other_contig, want=len(want), buckets={k: len(v) for k, v in buckets.items()})
json.dump(stats, open(f"{OUT}/{'smoke' if SMOKE else 'full'}.stats.json", "w"), indent=1)
log(json.dumps(stats))
print("MAP_SMOKE_OK" if SMOKE else "MAP_COMPLETE", flush=True)
