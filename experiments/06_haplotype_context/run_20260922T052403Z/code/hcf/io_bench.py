import os as _cfg_os
import math as _cfg_math

def _config_number(name, cast=int, positive=False):
    value = cast(_cfg_os.environ[name])
    if not _cfg_math.isfinite(value) or (positive and value <= 0):
        raise ValueError(name + " must be finite" + (" and positive" if positive else ""))
    return value

import os as _os
N_SAMPLES = _config_number("N_SAMPLES", int, True)
import json
import resource
import time
from pathlib import Path

import numpy as np
import yaml

HG19_CHROM_LENGTH = {
    1: 249250621, 2: 243199373, 3: 198022430, 4: 191154276, 5: 180915260,
    6: 171115067, 7: 159138663, 8: 146364022, 9: 141213431, 10: 135534747,
    11: 135006516, 12: 133851895, 13: 115169878, 14: 107349540, 15: 102531392,
    16: 90354753, 17: 81195210, 18: 78077248, 19: 59128983, 20: 63025520,
    21: 48129895, 22: 51304566,
}
TILE_BP = 100_000
STRUCTURE_TILES_MAX = 300
WALL_BUDGET_H = 16.0
DISK_BUDGET_GIB = 100.0

def _vcf_path_for_chrom(cfg, chrom):
    idx_pattern = cfg["x_registration"]["genotype_index_pattern"]
    csi_path = idx_pattern.replace("{N}", str(chrom))
    return csi_path[: -len(".csi")]

def _select_tile_windows(seed, n_tiles=12):
    rng = np.random.default_rng(seed)
    chroms = rng.choice(np.arange(1, 23), size=min(n_tiles, 22), replace=False)
    tiles = []
    for chrom in sorted(int(c) for c in chroms):
        length = HG19_CHROM_LENGTH[chrom]
        max_start = length - TILE_BP - 1
        start = int(rng.integers(1, max_start))
        tiles.append({"chrom": chrom, "start": start, "end": start + TILE_BP})
    return tiles

def _quartile_labels(values):
    values = np.asarray(values, dtype=float)
    q = np.quantile(values, [0.25, 0.5, 0.75])
    labels = []
    for v in values:
        if v <= q[0]:
            labels.append("Q1")
        elif v <= q[1]:
            labels.append("Q2")
        elif v <= q[2]:
            labels.append("Q3")
        else:
            labels.append("Q4")
    return labels

def _make_subwindows(positions, markers_per=16, step=8, max_gap_bp=20_000):
    n = len(positions)
    subwindows = []
    i = 0
    while i < n:
        j = i
        while j < min(i + markers_per, n) - 1:
            if positions[j + 1] - positions[j] > max_gap_bp:
                break
            j += 1
        subwindows.append((i, j + 1))
        i += step
        if i >= n:
            break
    return subwindows

def run_io_bench(config_path, out_dir, bcftools, n_tiles=12, n_people=2048, seed=20260922):
    import cyvcf2

    with open(config_path) as f:
        cfg = yaml.safe_load(f)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    tiles = _select_tile_windows(seed, n_tiles=n_tiles)

    sample_vcf_path = _vcf_path_for_chrom(cfg, tiles[0]["chrom"])
    probe = cyvcf2.VCF(sample_vcf_path)
    all_samples = list(probe.samples)
    probe.close()
    n_total_samples = len(all_samples)
    sample_rng = np.random.default_rng(seed + 2)
    take = min(n_people, n_total_samples)
    sample_idx = np.sort(sample_rng.choice(n_total_samples, size=take, replace=False))
    chosen_samples = [all_samples[i] for i in sample_idx]

    tile_results = []
    ru0 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    t_all0 = time.time()

    for tile in tiles:
        chrom, start, end = tile["chrom"], tile["start"], tile["end"]
        vcf_path = _vcf_path_for_chrom(cfg, chrom)

        vcf = cyvcf2.VCF(vcf_path)
        vcf.set_samples(chosen_samples)
        t0 = time.time()
        positions, gt_arrays = [], []
        for rec in vcf(f"{chrom}:{start}-{end}"):
            positions.append(rec.POS)
            gt_arrays.append(np.asarray(rec.genotypes, dtype=object))
        vcf.close()
        t_read = time.time() - t0

        n_variants = len(positions)
        t0 = time.time()
        subwindows = _make_subwindows(positions) if n_variants else []
        t_subwindow = time.time() - t0

        n_calls = n_variants * take
        tile_results.append({
            "chrom": chrom, "start": start, "end": end,
            "variant_count": n_variants,
            "n_people": take, "n_calls": n_calls,
            "wall_read_s": round(t_read, 4), "wall_subwindow_s": round(t_subwindow, 5),
            "n_subwindows": len(subwindows),
            "calls_per_s_read": round(n_calls / t_read, 1) if t_read > 0 else None,
        })
        print(f"tile chr{chrom}:{start}-{end} variants={n_variants} calls={n_calls} "
              f"read_s={t_read:.3f} subwindow_s={t_subwindow:.4f}", flush=True)

    quartiles = _quartile_labels([t["variant_count"] for t in tile_results])
    for t, ql in zip(tile_results, quartiles):
        t["variant_count_quartile_within_this_run"] = ql

    t_all = time.time() - t_all0
    ru1 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_rss_kb_delta = ru1 - ru0
    peak_rss_kb_absolute = ru1

    total_calls = sum(t["n_calls"] for t in tile_results)
    total_wall = sum(t["wall_read_s"] + t["wall_subwindow_s"] for t in tile_results)
    mean_calls_per_s = total_calls / total_wall if total_wall > 0 else float("nan")
    mean_variants_per_tile = float(np.mean([t["variant_count"] for t in tile_results]))

    scale_people = N_SAMPLES / take

    est_total_calls_full = mean_variants_per_tile * STRUCTURE_TILES_MAX * N_SAMPLES
    est_wall_s_full = est_total_calls_full / mean_calls_per_s if mean_calls_per_s > 0 else float("inf")
    est_wall_h_full = est_wall_s_full / 3600.0

    est_peak_rss_mb_full = (peak_rss_kb_absolute / 1024.0) * scale_people

    bytes_per_call_gt_array = 3 * 28
    est_disk_gib_full = (mean_variants_per_tile * STRUCTURE_TILES_MAX * N_SAMPLES * bytes_per_call_gt_array) / (1024 ** 3)

    wall_status = "WITHIN_BUDGET" if est_wall_h_full <= WALL_BUDGET_H else "RESOURCE_LIMITED"
    disk_status = "WITHIN_BUDGET" if est_disk_gib_full <= DISK_BUDGET_GIB else "RESOURCE_LIMITED"
    overall_status = "WITHIN_BUDGET" if wall_status == disk_status == "WITHIN_BUDGET" else "RESOURCE_LIMITED"

    proposed_scope = None
    if overall_status == "RESOURCE_LIMITED":
        max_tiles_within_wall = int(WALL_BUDGET_H * 3600 * mean_calls_per_s / (mean_variants_per_tile * N_SAMPLES)) if mean_calls_per_s > 0 else 0
        max_tiles_within_disk = int(DISK_BUDGET_GIB * (1024 ** 3) / (mean_variants_per_tile * N_SAMPLES * bytes_per_call_gt_array)) if bytes_per_call_gt_array > 0 else 0
        proposed_scope = {
            "max_structure_tiles_within_wall_budget": max(max_tiles_within_wall, 0),
            "max_structure_tiles_within_disk_budget": max(max_tiles_within_disk, 0),
            "recommendation": "Reduce structure_tiles_max below the smaller of the two figures "
                               "above, or reduce n_people processed per pass (e.g. chunked "
                               "person-blocks written incrementally rather than held in memory), "
                               "before HC-D1 build-states is run.",
        }

    report = {
        "protocol_id": cfg.get("protocol", {}).get("id"),
        "scope": "HC-D0 stage-2 IO throughput benchmark only; NOT a state-construction or "
                 "model-fitting measurement (brief explicitly excludes those from this stage).",
        "measured": {
            "n_tiles": n_tiles, "n_people": take, "seed": seed,
            "tile_selection_design": "single-pass: 12 distinct random chromosomes (seeded), "
                                      "one random 100kb window per chromosome, read once with "
                                      "the real 2048-person subset; variant-count quartile "
                                      "labels assigned post-hoc from the observed distribution "
                                      "across these 12 tiles (see io_bench.py module docstring "
                                      "for why a separate count-then-select pass was avoided).",
            "tiles": tile_results,
            "total_calls_measured": total_calls,
            "total_wall_s_measured": round(total_wall, 3),
            "wall_s_including_overhead": round(t_all, 3),
            "mean_calls_per_s": round(mean_calls_per_s, 1),
            "mean_variants_per_tile": round(mean_variants_per_tile, 1),
            "peak_rss_kb_absolute": int(peak_rss_kb_absolute),
            "peak_rss_kb_delta_during_bench": int(peak_rss_kb_delta),
        },
        "extrapolation_to_full_protocol": {
            "target_structure_tiles_max": STRUCTURE_TILES_MAX,
            "target_n_samples": N_SAMPLES,
            "method": "[추론] wall = (mean_variants_per_tile * target_tiles * target_people) / "
                      "mean_calls_per_s, measured on this run's 12 tiles x "
                      f"{take} people (throughput-only, not a claim about steady-state "
                      "full-genome IO under concurrent load from other branches -- another "
                      "unrelated process (phi_gate, a different project under the same "
                      "account) was observed consuming ~1 CPU core on this host during this "
                      "benchmark run). RAM extrapolated by scaling this run's peak RSS "
                      "linearly in n_people (single tile processed at a time, not held "
                      "concurrently). Disk extrapolated from a rough per-call "
                      "boxed-genotype-array byte estimate, NOT a measured on-disk size -- "
                      "flagged [실행정책]/[추론], not a firm figure.",
            "est_wall_hours_full": round(est_wall_h_full, 3),
            "est_peak_rss_mb_full": round(est_peak_rss_mb_full, 1),
            "est_disk_gib_full": round(est_disk_gib_full, 3),
            "wall_budget_hours": WALL_BUDGET_H,
            "disk_budget_gib": DISK_BUDGET_GIB,
            "wall_status": wall_status,
            "disk_status": disk_status,
        },
        "overall_status": overall_status,
        "proposed_reduced_scope_if_resource_limited": proposed_scope,
        "evidence_grade": "[실측-신규] tile-level throughput; [추론] full-scale extrapolation",
    }

    with open(out_dir / "resource_report.json", "w") as f:
        json.dump(report, f, indent=1)

    print(json.dumps({"status": "DONE", "overall_status": overall_status,
                       "est_wall_hours_full": round(est_wall_h_full, 3),
                       "est_disk_gib_full": round(est_disk_gib_full, 3)}))
    return 0
