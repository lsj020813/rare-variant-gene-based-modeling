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


import math, os, shutil, subprocess, sys, tempfile, glob

ROOT = _config_path("${PROJECT_ROOT}")
BCF = "bcftools"
PRODUCER = f"{ROOT}/scripts/measure_r2_af_fast_20260726.sh"
AGGREGATOR = f"{ROOT}/scripts/aggregate_g2_fmi_20260726.py"

FLOORS = ["0.9", "0.8", "0.5", "0.4", "0.3"]
RARE_MAX = 0.001

RECORDS = [
    (0.0002, 0.95), (0.0004, 0.85), (0.0006, 0.55), (0.0008, 0.45), (0.0010, 0.35),
    (0.0200, 0.99), (0.9995, 0.92), (0.3000, 0.25), (0.1000, 0.80), (0.0050, 0.30),
]

SPLIT = {
    "1": [(0.40, 0.98), (0.35, 0.97), (0.30, 0.99), (0.45, 0.96),
          (0.25, 0.98), (0.20, 0.97), (0.15, 0.99), (0.10, 0.98)],
    "2": [(0.0003, 0.32), (0.0004, 0.31)],
}

def m_of(af):
    return af if af <= 0.5 else 1 - af

def hand(records, floor, rare_only=False):
    kept = maf = info = 0
    for af, r2 in records:
        m = m_of(af)
        if rare_only and not (m < RARE_MAX):
            continue
        if r2 >= float(floor):
            kept += 1
            maf += m
            info += r2 * m
    if kept == 0:
        return kept, None, None, None, None
    frac = info / maf
    fmi = 1 - frac
    mraw = math.ceil(100 * fmi - 1e-12)
    return kept, maf, info, fmi, max(10, mraw)

def write_vcf(path, chrom, records):
    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write(f"##contig=<ID={chrom},length=250000000>\n")
        fh.write('##INFO=<ID=AF,Number=A,Type=Float,Description="af">\n')
        fh.write('##INFO=<ID=R2,Number=1,Type=Float,Description="r2">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for i, (af, r2) in enumerate(records, 1):
            fh.write(f"{chrom}\t{1000*i}\t.\tA\tG\t.\tPASS\tAF={af:.17g};R2={r2:.17g}\n")
    subprocess.run([BCF, "view", "-Oz", "-o", path + ".gz", path], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    subprocess.run([BCF, "index", "-t", path + ".gz"], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return path + ".gz"

def read_kv(path):
    kv = {}
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) == 2:
                kv[p[0]] = p[1]
    return kv

def gate_producer(work):
    print("=" * 78)
    print("GATE 1 — producer (measure_r2_af_fast): 10 레코드 손계산 vs 스크립트")
    print("=" * 78)
    vcf = write_vcf(os.path.join(work, "g1.vcf"), "21", RECORDS)
    out = os.path.join(work, "g1_out.tsv")
    tbl = os.path.join(work, "g1_table.tsv.gz")
    r = subprocess.run(["bash", PRODUCER, "--chrom", "21", "--vcf", vcf,
                        "--out", out, "--emit-table", tbl],
                       capture_output=True, text=True)
    if r.returncode != 0:
        print("  스크립트 실행 실패:", (r.stderr or r.stdout).strip()[:400])
        return False
    kv = read_kv(out)

    ok = True
    for band, pref, rare in (("전체", "floor_", False), ("희귀 m<0.001", "rare_floor_", True)):
        print(f"\n  [{band}]")
        print(f"  {'floor':>5} {'kept':>5} {'손계산 FMI':>22} {'M':>3} | {'스크립트 FMI':>22} {'M':>3} | 차이     일치")
        for f in FLOORS:
            k, ms, isum, fmi, M = hand(RECORDS, f, rare)
            gk = int(kv[f"{pref}{f}_kept"])
            if k == 0:
                good = gk == 0
                print(f"  {f:>5} {k:>5} {'(빈 층)':>22} {'-':>3} | {'':>22} {'-':>3} | {'':>8} {'O' if good else 'X'}")
                ok &= good
                continue
            gfmi = float(kv[f"{pref}{f}_FMI"])
            gM = int(kv[f"{pref}{f}_White_M_min10"])
            gms = float(kv[f"{pref}{f}_maf_sum"])
            gis = float(kv[f"{pref}{f}_info_sum"])
            d = abs(fmi - gfmi)
            good = (k == gk and M == gM and d == 0.0
                    and abs(ms - gms) == 0.0 and abs(isum - gis) == 0.0)
            ok &= good
            print(f"  {f:>5} {k:>5} {fmi:>22.17f} {M:>3} | {gfmi:>22.17f} {gM:>3} | {d:.1e} {'O' if good else 'X'}")

    import gzip
    rows = [l.split("\t") for l in gzip.open(tbl, "rt").read().splitlines()[1:]]
    tb = [(float(r[4]), float(r[5])) for r in rows]
    for f in FLOORS:
        _, _, _, fmi, _ = hand(tb, f)
        good = fmi == float(kv[f"floor_{f}_FMI"])
        ok &= good
        if not good:
            print(f"  emit-table 재계산 floor {f} 불일치")
    print(f"\n  emit-table 재계산으로 전 floor 요약 재현: {'O' if ok else 'X'}")
    return ok

def gate_aggregator(work):
    print()
    print("=" * 78)
    print("GATE 2 — aggregator: 같은 10 레코드를 2 염색체에 8:2 분할")
    print("=" * 78)
    rundir = os.path.join(work, "aggrun")
    os.makedirs(os.path.join(rundir, "out"), exist_ok=True)
    for c, recs in SPLIT.items():
        vcf = write_vcf(os.path.join(work, f"g2_c{c}.vcf"), c, recs)
        r = subprocess.run(["bash", PRODUCER, "--chrom", c, "--vcf", vcf,
                            "--out", os.path.join(rundir, "out", f"chr{c}.tsv")],
                           capture_output=True, text=True)
        if r.returncode != 0:
            print("  producer 실패:", (r.stderr or r.stdout).strip()[:300])
            return False

    import importlib.util
    spec = importlib.util.spec_from_file_location("aggmod", AGGREGATOR)
    agg = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(agg)
    agg.EXPECTED_CHROMS = ["1", "2"]
    agg.EXPECTED_TOTAL_RECORDS = sum(len(v) for v in SPLIT.values())

    import io, contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        agg.main(rundir)
    out = buf.getvalue()

    allrec = [r for v in SPLIT.values() for r in v]
    ok = True
    print(f"  {'floor':>5} {'손계산 FMI':>22} {'M':>3} | {'집계기 FMI':>22} {'M':>3} | 차이     일치")
    for f in FLOORS:
        k, ms, isum, fmi, M = hand(allrec, f)
        row = [l for l in out.splitlines() if l.startswith(f + "\t")]
        if not row:
            print(f"  {f:>5} 집계기 출력 없음"); ok = False; continue
        cols = row[0].split("\t")
        gfmi, gM = float(cols[5]), int(cols[7])
        d = abs(fmi - gfmi)
        good = (d == 0.0 and M == gM)
        ok &= good
        print(f"  {f:>5} {fmi:>22.17f} {M:>3} | {gfmi:>22.17f} {gM:>3} | {d:.1e} {'O' if good else 'X'}")

    f = "0.3"
    per = [hand(v, f)[3] for v in SPLIT.values()]
    wrong = sum(per) / len(per)
    right = hand(allrec, f)[3]
    cols = [l for l in out.splitlines() if l.startswith(f + "\t")][0].split("\t")
    got = float(cols[5])
    sep = abs(right - wrong)
    print(f"\n  금지된 방식(염색체 FMI 평균) = {wrong:.8f} → M={max(10, math.ceil(100*wrong-1e-12))}")
    print(f"  올바른 방식(Σinfo/Σmaf)      = {right:.8f} → M={max(10, math.ceil(100*right-1e-12))}")
    print(f"  두 방식 간격 {sep*100:.2f}%p — 시험이 둘을 구분할 수 있는가: {'O' if sep > 0.01 else 'X 시험 무의미'}")
    disc = abs(got - wrong) > 1e-6
    ok &= disc and sep > 0.01
    print(f"  집계기가 금지된 방식이 아님: {'O' if disc else 'X'}")
    return ok

class Tee:
    def __init__(self, *streams): self.streams = streams
    def write(self, s):
        for st in self.streams: st.write(s)
    def flush(self):
        for st in self.streams: st.flush()

def sha(path):
    import hashlib
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()

def main():
    import datetime, io, contextlib
    keep = "--keep" in sys.argv
    logdir = f"{ROOT}/reports/gate10"
    os.makedirs(logdir, exist_ok=True)
    ts = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    logpath = f"{logdir}/gate10_{ts}.log"
    work = tempfile.mkdtemp(prefix="gate10_", dir=_config_path("${PROJECT_ROOT}/work/tmp"))
    rc = 1
    with open(logpath, "w") as lf:
        tee = Tee(sys.stdout, lf)
        with contextlib.redirect_stdout(tee):
            print(f"gate10 run          : {ts}")
            print(f"gate script sha256  : {sha(__file__)}")
            print(f"producer sha256     : {sha(PRODUCER)}")
            print(f"aggregator sha256   : {sha(AGGREGATOR)}")
            print(f"bcftools            : {subprocess.run([BCF,'--version'],capture_output=True,text=True).stdout.splitlines()[0]}")
            print(f"records under test  : {len(RECORDS)}")
            for i, (af, r2) in enumerate(RECORDS, 1):
                print(f"  rec{i:>2}  AF={af:<10g} R2={r2:<6g} m=min(AF,1-AF)={m_of(af):<10g} rare={'yes' if m_of(af) < RARE_MAX else 'no'}")
            print()
            try:
                r1 = gate_producer(work)
                r2 = gate_aggregator(work)
                print()
                print("=" * 78)
                print(f"GATE 1 producer   : {'PASS' if r1 else 'FAIL'}")
                print(f"GATE 2 aggregator : {'PASS' if r2 else 'FAIL'}")
                print(f"판정: {'스크립트 완결 — 10개 데이터 손계산과 정확히 일치' if (r1 and r2) else '미완결 — 불일치 항목 위 참조'}")
                print("=" * 78)
                rc = 0 if (r1 and r2) else 1
            finally:
                if keep:
                    print(f"scratch: {work}")
                else:
                    shutil.rmtree(work, ignore_errors=True)
    print(f"\n로그 기록: {logpath}")
    return rc

if __name__ == "__main__":
    sys.exit(main())
