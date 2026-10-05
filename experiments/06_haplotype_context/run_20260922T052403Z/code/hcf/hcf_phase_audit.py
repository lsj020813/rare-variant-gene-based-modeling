#!/usr/bin/env python
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

import os as _os
N_SAMPLES = _config_number("N_SAMPLES", int, True)
import subprocess, sys, os, json, time, math, hashlib, csv
import numpy as np
from scipy.stats import binomtest

BCF = "bcftools"
OI = _config_path("${PROJECT_ROOT}/work/ref/orig_index")
S37 = _config_path("${PROJECT_ROOT}/work/tmp/lift38/extracts")
OUT = sys.argv[1]
CHRS = [int(x) for x in sys.argv[2].split(",")]
MODE = sys.argv[3] if len(sys.argv) > 3 else "audit"
TAG = sys.argv[4] if len(sys.argv) > 4 else ""
SEED = 20260922
NREC = 2000
NSAMP = 1000
PROTO = "HCF-20260922-v1"
NICE = ["nice", "-n", "19", "ionice", "-c3"]
os.makedirs(f"{OUT}/tmp", exist_ok=True)

def log(msg):
    with open(f"{OUT}/phase_audit{TAG}.log", "a") as f:
        f.write(time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()) + " " + msg + "\n")

def bq(args):
    p = subprocess.run(NICE + [BCF] + args, capture_output=True)
    if p.returncode != 0:
        raise RuntimeError(p.stderr.decode()[-2000:])
    return p.stdout.decode()

def load_sites(chrN):
    txt = bq(["query", "-f", "%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/AVG_CS\t%INFO/R2\t%INFO/TYPED\n", f"{S37}/chr{chrN}.s37.vcf.gz"])
    lines = txt.rstrip("\n").split("\n")
    return lines

def sample_ids_select():
    ids = bq(["query", "-l", f"{OI}/chr22.vcf.gz"]).split()
    assert len(ids) == N_SAMPLES, len(ids)
    rng = np.random.RandomState(SEED)
    sidx = np.sort(rng.choice(len(ids), NSAMP, replace=False))
    sel = [ids[i] for i in sidx]
    md5 = hashlib.md5(("\n".join(sel)).encode()).hexdigest()
    return sel, md5

def audit_chr(chrN, sample_ids, smd5):
    t0 = time.time()
    lines = load_sites(chrN)
    n_sites = len(lines)
    rng = np.random.RandomState(SEED * 100 + chrN)
    ridx = np.sort(rng.choice(n_sites, NREC, replace=False))
    chosen = {}
    positions = set()
    for i in ridx:
        pos, ref, alt, af, cs, r2, typed = lines[i].split("\t")
        chosen[(pos, ref, alt)] = (af, cs, r2)
        positions.add(int(pos))
    reg = f"{OUT}/tmp/regions_chr{chrN}.tsv"
    with open(reg, "w") as f:
        for p in sorted(positions):
            f.write(f"{chrN}\t{p}\n")
    fmt = "%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/AVG_CS\t%INFO/R2\t%INFO/TYPED[\t%GT:%GP:%DS]\n"
    cmd = NICE + [BCF, "query", "-R", reg, "-S", "/dev/stdin", "-f", fmt, f"{OI}/chr{chrN}.vcf.gz"]
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc.stdin.write(("\n".join(sample_ids) + "\n").encode())
    proc.stdin.close()
    st = dict(n_records_returned=0, n_records_matched=0, n_extra_records_same_pos=0, info_match_s37=0,
              n_calls=0, n_phased=0, n_unphased=0, n_missing=0, n_haploid_or_other=0,
              agree_ds=0, agree_gp=0, gp_ties=0, ds_eq_gp=0,
              typed_records=0, imputed_records=0,
              calls_typed=0, calls_imputed=0, agree_ds_typed=0, agree_ds_imputed=0, agree_gp_typed=0, agree_gp_imputed=0,
              n_het=0, n_01=0, n_10=0, n_het_other=0,
              n_records_het_ge20=0, n_records_direction_imbalanced_p001=0,
              n_gt_1_gp_disagree_at_het=0)
    ps01 = np.zeros(NSAMP, dtype=np.int64)
    pshet = np.zeros(NSAMP, dtype=np.int64)
    ns_seen = None
    for raw in proc.stdout:
        line = raw.decode().rstrip("\n")
        parts = line.split("\t")
        pos, ref, alt, af, cs, r2, typed = parts[:7]
        calls = parts[7:]
        st["n_records_returned"] += 1
        key = (pos, ref, alt)
        if key not in chosen:
            st["n_extra_records_same_pos"] += 1
            continue
        st["n_records_matched"] += 1
        if chosen[key] == (af, cs, r2):
            st["info_match_s37"] += 1
        is_typed = (typed == "1")
        if is_typed:
            st["typed_records"] += 1
        else:
            st["imputed_records"] += 1
        if ns_seen is None:
            ns_seen = len(calls)
        c01 = c10 = 0
        for j, cstr in enumerate(calls):
            gt, gp, ds = cstr.split(":")
            st["n_calls"] += 1
            if "." in gt:
                st["n_missing"] += 1
                continue
            if "|" in gt:
                st["n_phased"] += 1
                a, b = gt.split("|")
            elif "/" in gt:
                st["n_unphased"] += 1
                a, b = gt.split("/")
            else:
                st["n_haploid_or_other"] += 1
                continue
            a = int(a); b = int(b); s = a + b
            d = float(ds)
            rd = int(math.floor(d + 0.5))
            gpf = [float(x) for x in gp.split(",")]
            mx = max(gpf)
            am = gpf.index(mx)
            if gpf.count(mx) > 1:
                st["gp_ties"] += 1
            ok_ds = (s == rd); ok_gp = (am == s)
            if ok_ds: st["agree_ds"] += 1
            if ok_gp: st["agree_gp"] += 1
            if abs(d - (gpf[1] + 2 * gpf[2])) <= 0.0105:
                st["ds_eq_gp"] += 1
            if is_typed:
                st["calls_typed"] += 1
                if ok_ds: st["agree_ds_typed"] += 1
                if ok_gp: st["agree_gp_typed"] += 1
            else:
                st["calls_imputed"] += 1
                if ok_ds: st["agree_ds_imputed"] += 1
                if ok_gp: st["agree_gp_imputed"] += 1
            if s == 1 and "|" in gt:
                st["n_het"] += 1
                pshet[j] += 1
                if not ok_gp:
                    st["n_gt_1_gp_disagree_at_het"] += 1
                if gt == "0|1":
                    st["n_01"] += 1; c01 += 1; ps01[j] += 1
                elif gt == "1|0":
                    st["n_10"] += 1; c10 += 1
                else:
                    st["n_het_other"] += 1
        if c01 + c10 >= 20:
            st["n_records_het_ge20"] += 1
            if binomtest(c01, c01 + c10, 0.5).pvalue < 0.001:
                st["n_records_direction_imbalanced_p001"] += 1
    err = proc.stderr.read().decode()
    rc = proc.wait()
    if rc != 0:
        raise RuntimeError(f"bcftools query rc={rc}: {err[-2000:]}")
    mask = pshet >= 20
    p01 = ps01[mask] / pshet[mask]
    n_biased = 0
    for k in np.where(mask)[0]:
        if binomtest(int(ps01[k]), int(pshet[k]), 0.5).pvalue < 0.001:
            n_biased += 1
    nc = st["n_calls"]; nv = st["n_phased"] + st["n_unphased"]
    het = st["n_01"] + st["n_10"]
    row = dict(
        protocol_id=PROTO, chromosome=chrN,
        source_path=f"{OI}/chr{chrN}.vcf.gz", realpath=os.path.realpath(f"{OI}/chr{chrN}.vcf.gz"),
        sampling_frame=f"{S37}/chr{chrN}.s37.vcf.gz (bcftools view -G of realpath, sites-only)",
        n_records_frame=n_sites,
        sampling_rule="uniform random records without replacement from sites frame; uniform random samples without replacement from header sample list (same 1000 for all chromosomes)",
        seed_samples=SEED, seed_records=SEED * 100 + chrN, rng="numpy.RandomState",
        n_records_requested=NREC, n_samples_requested=NSAMP, n_samples_returned=ns_seen if ns_seen else 0,
        sample_selection_md5=smd5,
        n_records_returned=st["n_records_returned"], n_records_matched=st["n_records_matched"],
        n_extra_records_same_pos=st["n_extra_records_same_pos"],
        info_match_s37_frac=round(st["info_match_s37"] / max(1, st["n_records_matched"]), 6),
        typed_records=st["typed_records"], imputed_records=st["imputed_records"],
        n_calls=nc,
        phased_frac=round(st["n_phased"] / max(1, nc), 6),
        unphased_frac=round(st["n_unphased"] / max(1, nc), 6),
        missing_frac=round(st["n_missing"] / max(1, nc), 6),
        other_gt_frac=round(st["n_haploid_or_other"] / max(1, nc), 6),
        gt_eq_round_ds_frac=round(st["agree_ds"] / max(1, nv), 6),
        gt_eq_argmax_gp_frac=round(st["agree_gp"] / max(1, nv), 6),
        gp_tie_frac=round(st["gp_ties"] / max(1, nv), 6),
        ds_eq_gp1_plus_2gp2_frac=round(st["ds_eq_gp"] / max(1, nv), 6),
        gt_eq_round_ds_frac_typed=round(st["agree_ds_typed"] / max(1, st["calls_typed"]), 6) if st["calls_typed"] else "",
        gt_eq_round_ds_frac_imputed=round(st["agree_ds_imputed"] / max(1, st["calls_imputed"]), 6) if st["calls_imputed"] else "",
        gt_eq_argmax_gp_frac_typed=round(st["agree_gp_typed"] / max(1, st["calls_typed"]), 6) if st["calls_typed"] else "",
        gt_eq_argmax_gp_frac_imputed=round(st["agree_gp_imputed"] / max(1, st["calls_imputed"]), 6) if st["calls_imputed"] else "",
        n_het=st["n_het"], n_01=st["n_01"], n_10=st["n_10"], n_het_other=st["n_het_other"],
        frac_01=round(st["n_01"] / max(1, het), 6),
        binom_p_01_vs_10=(binomtest(st["n_01"], het, 0.5).pvalue if het else ""),
        het_gp_argmax_disagree_frac=round(st["n_gt_1_gp_disagree_at_het"] / max(1, st["n_het"]), 6),
        per_sample_p01_mean=round(float(p01.mean()), 6) if len(p01) else "",
        per_sample_p01_sd=round(float(p01.std()), 6) if len(p01) else "",
        per_sample_p01_min=round(float(p01.min()), 6) if len(p01) else "",
        per_sample_p01_max=round(float(p01.max()), 6) if len(p01) else "",
        n_samples_het_ge20=int(mask.sum()), n_samples_direction_biased_p001=n_biased,
        n_records_het_ge20=st["n_records_het_ge20"],
        n_records_direction_imbalanced_p001=st["n_records_direction_imbalanced_p001"],
        wall_s=round(time.time() - t0, 1),
    )
    log(f"chr{chrN} done: " + json.dumps({k: row[k] for k in ['n_records_matched','n_calls','phased_frac','missing_frac','gt_eq_round_ds_frac','gt_eq_argmax_gp_frac','frac_01','wall_s']}))
    return row

def chunk_scan(chrN):
    t0 = time.time()
    lines = load_sites(chrN)
    n = len(lines)
    pos = np.empty(n, dtype=np.int64); cs = np.full(n, np.nan); r2 = np.full(n, np.nan); typed = np.zeros(n, dtype=bool)
    keys = []
    for i, l in enumerate(lines):
        p, ref, alt, af, c, r, t = l.split("\t")
        pos[i] = int(p)
        try: cs[i] = float(c)
        except ValueError: pass
        try: r2[i] = float(r)
        except ValueError: pass
        typed[i] = (t == "1")
        keys.append((p, ref, alt))
    dup_keys = n - len(set(keys))
    sorted_ok = bool(np.all(np.diff(pos) >= 0))
    d = np.diff(pos)
    gap_thr = 500_000
    gi = np.where(d > gap_thr)[0]
    gaps = [dict(pos_before=int(pos[i]), pos_after=int(pos[i + 1]), gap_bp=int(d[i])) for i in gi]
    W = 1000; step = 250
    drops = []
    cands = []
    for b in range(W, n - W, step):
        mb_r = np.nanmedian(r2[b - W:b]); ma_r = np.nanmedian(r2[b:b + W])
        mb_c = np.nanmedian(cs[b - W:b]); ma_c = np.nanmedian(cs[b:b + W])
        dr = ma_r - mb_r; dc = ma_c - mb_c
        if abs(dr) > 0.15 or abs(dc) > 0.03:
            cands.append(dict(pos=int(pos[b]), record_index=int(b), r2_med_before=round(float(mb_r), 4), r2_med_after=round(float(ma_r), 4),
                              cs_med_before=round(float(mb_c), 4), cs_med_after=round(float(ma_c), 4)))
    merged = []
    for c in cands:
        if merged and c["record_index"] - merged[-1]["record_index_last"] <= 2000:
            merged[-1]["record_index_last"] = c["record_index"]; merged[-1]["pos_last"] = c["pos"]
            merged[-1]["max_abs_r2_step"] = max(merged[-1]["max_abs_r2_step"], abs(c["r2_med_after"] - c["r2_med_before"]))
            merged[-1]["max_abs_cs_step"] = max(merged[-1]["max_abs_cs_step"], abs(c["cs_med_after"] - c["cs_med_before"]))
        else:
            merged.append(dict(pos_first=c["pos"], pos_last=c["pos"], record_index_first=c["record_index"], record_index_last=c["record_index"],
                               max_abs_r2_step=abs(c["r2_med_after"] - c["r2_med_before"]), max_abs_cs_step=abs(c["cs_med_after"] - c["cs_med_before"])))
    hyp = []
    for m in range(20_000_000, int(pos.max()), 20_000_000):
        near = (np.abs(pos - m) <= 250_000); flank = (np.abs(pos - m) > 250_000) & (np.abs(pos - m) <= 1_000_000)
        if near.sum() >= 50 and flank.sum() >= 50:
            hyp.append(dict(boundary=m, n_near=int(near.sum()), n_flank=int(flank.sum()),
                            r2_mean_near=round(float(np.nanmean(r2[near])), 4), r2_mean_flank=round(float(np.nanmean(r2[flank])), 4),
                            cs_mean_near=round(float(np.nanmean(cs[near])), 4), cs_mean_flank=round(float(np.nanmean(cs[flank])), 4)))
        else:
            hyp.append(dict(boundary=m, n_near=int(near.sum()), n_flank=int(flank.sum()), note="insufficient records (gap/centromere)"))
    hdr = bq(["view", "-h", f"{OI}/chr{chrN}.vcf.gz"])
    hl = [l for l in hdr.split("\n") if l.startswith("##") and not l.startswith(("##INFO", "##FORMAT", "##contig"))]
    contig = [l for l in hdr.split("\n") if l.startswith("##contig")]
    clue = [l for l in hdr.split("\n") if l.startswith("##") and any(k in l.lower() for k in ["chunk", "overlap", "region", "command", "minimac4_", "pipeline", "window", "start", "end="])]
    res = dict(protocol_id=PROTO, chromosome=chrN, frame=f"{S37}/chr{chrN}.s37.vcf.gz", n_records=n, duplicate_keys=dup_keys, position_sorted=sorted_ok,
               pos_min=int(pos.min()), pos_max=int(pos.max()), n_typed=int(typed.sum()),
               r2_quantiles={q: round(float(np.nanquantile(r2, q)), 4) for q in [0.01, 0.1, 0.5, 0.9]},
               cs_quantiles={q: round(float(np.nanquantile(cs, q)), 4) for q in [0.01, 0.1, 0.5, 0.9]},
               gap_threshold_bp=gap_thr, n_gaps=len(gaps), gaps=gaps[:60],
               rolling_window_records=W, rolling_step=step, thresholds=dict(r2_median_step=0.15, cs_median_step=0.03),
               n_rolling_candidates=len(cands), n_candidate_loci_merged=len(merged), candidate_loci=merged[:60],
               minimac4_default_chunk_hypothesis_20Mb=hyp,
               header_non_info_lines=hl, header_contig_lines=contig, header_chunk_clue_lines=clue,
               wall_s=round(time.time() - t0, 1))
    with open(f"{OUT}/chunk_scan_chr{chrN}.json", "w") as f:
        json.dump(res, f, indent=1)
    log(f"chunkscan chr{chrN} done n={n} gaps={len(gaps)} cand={len(merged)} wall={res['wall_s']}")
    return res

if __name__ == "__main__":
    log(f"start mode={MODE} chrs={CHRS} pid={os.getpid()}")
    if MODE == "audit":
        sids, smd5 = sample_ids_select()
        log(f"samples selected n={len(sids)} md5={smd5}")
        outcsv = f"{OUT}/phase_readiness_raw{TAG}.csv"
        rows = []
        if os.path.exists(outcsv):
            with open(outcsv) as f:
                rows = list(csv.DictReader(f))
        done = {int(r["chromosome"]) for r in rows}
        for chrN in CHRS:
            if chrN in done:
                continue
            try:
                row = audit_chr(chrN, sids, smd5)
            except Exception as e:
                log(f"chr{chrN} ERROR {e}")
                row = dict(protocol_id=PROTO, chromosome=chrN, error=str(e)[:500])
            rows.append(row)
            keys = []
            for r in rows:
                for k in r:
                    if k not in keys: keys.append(k)
            with open(outcsv + ".tmp", "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)
            os.replace(outcsv + ".tmp", outcsv)
        open(f"{OUT}/audit{TAG}.done", "w").write(time.strftime("%Y%m%dT%H%M%SZ", time.gmtime()))
    else:
        for chrN in CHRS:
            try:
                chunk_scan(chrN)
            except Exception as e:
                log(f"chunkscan chr{chrN} ERROR {e}")
        open(f"{OUT}/chunkscan.done", "w").write(time.strftime("%Y%m%dT%H%M%SZ", time.gmtime()))
    log("end")
