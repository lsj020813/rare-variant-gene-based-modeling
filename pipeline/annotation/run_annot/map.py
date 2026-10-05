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


import gzip, glob, os, subprocess, sys, collections

R    = _config_path("${PROJECT_ROOT}/work/ref")
GP   = f"{R}/groupfiles_bwg"
KEY  = f"{R}/lift38_keyed"
OUT  = f"{R}/annot/map38"
BCF  = "bcftools"
SMOKE = "--smoke" in sys.argv
os.makedirs(OUT, exist_ok=True)

CHRS = [22] if SMOKE else list(range(1, 23))

want = set()
for n in range(1, 23):
    fp = f"{GP}/chr{n}.B_3kb_re2g.txt"
    if not os.path.exists(fp): continue
    with open(fp) as fh:
        for line in fh:
            f = line.rstrip("\n").split()
            if len(f) > 2 and f[1] == "var":
                want.update(f[2:])
print(f"[map] 그룹파일 참조 37키: {len(want):,}", flush=True)
assert want, "GATE FAIL: 그룹파일에서 키를 하나도 읽지 못했다"

buckets = collections.defaultdict(list)
seen37 = set()
src = CHRS if SMOKE else range(1, 23)
for n in src:
    fp = f"{KEY}/chr{n}.keyed38.vcf.gz"
    if not os.path.exists(fp):
        print(f"[map] 없음: {fp}", flush=True); continue
    p = subprocess.Popen(
        [BCF, "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t%ID\n", fp],
        stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    kept = tot = 0
    for line in p.stdout:
        tot += 1
        ch, pos, ref, alt, vid = line.rstrip("\n").split("\t")
        k37 = vid[3:] if vid.startswith("chr") else vid
        if k37 not in want: continue
        cb = ch[3:] if ch.startswith("chr") else ch
        buckets[cb].append((f"{cb}:{pos}:{ref}:{alt}", k37))
        seen37.add(k37); kept += 1
        if SMOKE and kept >= 200000: break
    p.stdout.close(); p.wait()
    print(f"[map] chr{n}: 주행 {tot:,}  채택 {kept:,}", flush=True)

assert seen37, "GATE FAIL: 38↔37 대응이 0건 — 키 형식 불일치"
if not SMOKE:
    cov = len(seen37) / len(want)
    print(f"[map] 대응 확보 37키 {len(seen37):,} / {len(want):,} = {cov*100:.2f}%", flush=True)
    assert cov > 0.90, f"GATE FAIL: 대응률 {cov*100:.1f}% — 90% 미만"

for cb, rows in sorted(buckets.items(), key=lambda kv: (len(kv[0]), kv[0])):
    with open(f"{OUT}/chr{cb}.map.tsv", "w") as fh:
        for k38, k37 in rows:
            fh.write(f"{k38}\t{k37}\n")
    print(f"[map] chr{cb}.map.tsv  {len(rows):,}행", flush=True)
print("MAP_SMOKE_OK" if SMOKE else "MAP_COMPLETE", flush=True)
