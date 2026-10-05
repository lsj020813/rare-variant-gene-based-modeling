
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import argparse
import json
import os
import resource
import subprocess
import sys
import threading
import time

import numpy as np
import pandas as pd

from hcf import build_states as bs
from hcf import fastpath
from hcf import tiles as htiles
from hcf import vcfio

ORIG = _config_path("${PROJECT_ROOT}/work/ref/orig_index/chr%s.vcf.gz")
WALL_BUDGET_H = 16.0
DISK_BUDGET_GIB = 100.0

def _maxrss_gib():
    a = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    b = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    return max(a, b) / 1024.0 / 1024.0

class SamplePipe(object):

    def __init__(self, ids):
        self.ids = list(ids)

    def __enter__(self):
        self.r, self.w = os.pipe()
        payload = ("\n".join(self.ids) + "\n").encode("ascii")

        def writer():
            try:
                os.write(self.w, payload)
            finally:
                os.close(self.w)

        self.th = threading.Thread(target=writer)
        self.th.daemon = True
        self.th.start()
        return "/dev/fd/%d" % self.r, (self.r,)

    def __exit__(self, *a):
        self.th.join(timeout=5)
        try:
            os.close(self.r)
        except OSError:
            pass
        return False

def bench_tile(chrom, tile_start, vt, n_samples, sample_ids=None):
    keep = set(zip(vt["pos"].tolist(), vt["ref"].tolist(), vt["alt"].tolist()))
    region = "%s:%d-%d" % (chrom, tile_start + 1, tile_start + htiles.TILE_BP)
    t0 = time.time()
    if sample_ids is None:
        rd = vcfio.read_region(ORIG % chrom, region, n_samples=n_samples, keep=keep)
    else:
        with SamplePipe(sample_ids) as (spath, fds):
            rd = vcfio.read_region(ORIG % chrom, region, n_samples=n_samples,
                                   samples_file=spath, keep=keep, pass_fds=fds)
    t_read = time.time() - t0
    hap = vcfio.to_person_major(rd["hap"])
    t1 = time.time()
    sw = htiles.subwindows(rd["pos"], markers=16, step=8, max_gap_bp=20000)
    n_keys = 0
    for (a, b) in sw:
        k = fastpath.encode_keys(hap[:, :, a:b])
        n_keys += k.size
    t_sub = time.time() - t1
    return {
        "chrom": str(chrom), "tile_start": int(tile_start),
        "tile_id": "chr%s:%d" % (chrom, tile_start),
        "n_people": int(n_samples),
        "n_records_in_tile": int(rd["n_lines"]),
        "n_primary_variants": int(len(vt)),
        "n_variants_read": int(rd["hap"].shape[0]),
        "n_subwindows": len(sw),
        "read_wall_s": round(t_read, 3),
        "subwindow_build_wall_s": round(t_sub, 3),
        "total_wall_s": round(t_read + t_sub, 3),
        "genotype_calls": int(rd["hap"].shape[0]) * int(n_samples),
        "calls_per_s": (int(rd["hap"].shape[0]) * int(n_samples) / t_read) if t_read > 0 else float("nan"),
        "hap_array_mib": float(hap.nbytes) / 1048576.0,
        "sep_all_pipe": bool(rd["sep_all_pipe"]),
    }

def bench_view_form(chrom, tile_start, n_samples, sample_ids):
    region = "%s:%d-%d" % (chrom, tile_start + 1, tile_start + htiles.TILE_BP)
    with SamplePipe(sample_ids) as (spath, fds):
        t0 = time.time()
        p = subprocess.Popen([vcfio.BCFTOOLS, "view", "-H", "-r", region, "-S", spath,
                              ORIG % chrom], stdout=subprocess.PIPE, pass_fds=fds)
        nb = 0
        nl = 0
        for raw in p.stdout:
            nb += len(raw)
            nl += 1
        p.stdout.close()
        p.wait()
        t = time.time() - t0
    return {"tile_id": "chr%s:%d" % (chrom, tile_start), "n_people": int(n_samples),
            "view_wall_s": round(t, 3), "view_bytes": nb, "view_records": nl}

def main(argv=None):
    ap = argparse.ArgumentParser("hcf.bench_io")
    ap.add_argument("--bench-tiles", required=True)
    ap.add_argument("--vardir", required=True)
    ap.add_argument("--split", required=True)
    ap.add_argument("--sample-list", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tmpdir", required=True)
    ap.add_argument("--full-probe-tiles", type=int, default=3)
    ap.add_argument("--bench-people", type=int, default=2048)
    ap.add_argument("--structure-tiles", type=int, default=300)
    args = ap.parse_args(argv)

    t_start = time.time()
    ids = [l.strip() for l in open(args.sample_list) if l.strip()]
    n_all = len(ids)
    rs = np.random.RandomState(20260922)
    bench_ids = [ids[i] for i in sorted(rs.choice(n_all, size=args.bench_people, replace=False))]

    sp = pd.read_csv(args.split, sep="\t")
    splits = dict((k, sp.index[sp["split"] == k].values) for k in ("A", "B", "C"))
    batch_proxy = sp["batch_proxy"].values

    bt = pd.read_csv(args.bench_tiles, sep="\t", dtype={"chrom": str})
    vcache = {}
    rows_small, rows_full, rows_view, rows_pipeline = [], [], [], []
    for i in range(len(bt)):
        ch = str(bt["chrom"].iloc[i])
        ts = int(bt["tile_start"].iloc[i])
        if ch not in vcache:
            vcache[ch] = bs.load_variants(ch, args.vardir)
        va = vcache[ch]
        vt = va[(va["pos"] >= ts) & (va["pos"] < ts + htiles.TILE_BP)].reset_index(drop=True)
        r = bench_tile(ch, ts, vt, args.bench_people, sample_ids=bench_ids)
        r["quartile"] = int(bt["variant_count_quartile"].iloc[i])
        rows_small.append(r)
        sys.stderr.write("small %s %s\n" % (r["tile_id"], r["total_wall_s"]))
        sys.stderr.flush()

    for i in range(min(args.full_probe_tiles, len(bt))):
        ch = str(bt["chrom"].iloc[i])
        ts = int(bt["tile_start"].iloc[i])
        va = vcache[ch]
        vt = va[(va["pos"] >= ts) & (va["pos"] < ts + htiles.TILE_BP)].reset_index(drop=True)
        rows_full.append(bench_tile(ch, ts, vt, n_all, sample_ids=None))
        rows_view.append(bench_view_form(ch, ts, args.bench_people, bench_ids))
        out = {"state_inventory": [], "diplotype_contrast": [], "phase_pairs": [],
               "repr_info": [], "tile_status": [], "tmpdir": args.tmpdir}
        t0 = time.time()
        st = bs.process_tile(ch, ts, vt, splits, batch_proxy, out, verify=True,
                             n_people_expected=n_all)
        st["pipeline_wall_s"] = round(time.time() - t0, 2)
        st["rows_state_inventory"] = len(out["state_inventory"])
        st["rows_diplotype"] = len(out["diplotype_contrast"])
        st["rows_phase_pairs"] = len(out["phase_pairs"])
        st["rows_repr_info"] = len(out["repr_info"])
        rows_pipeline.append(st)
        sys.stderr.write("full %s pipeline %.1fs\n" % (st["tile_id"], st["pipeline_wall_s"]))
        sys.stderr.flush()

    small = pd.DataFrame(rows_small)
    full = pd.DataFrame(rows_full)
    pipe = pd.DataFrame(rows_pipeline)

    per_tile_pipeline = float(pipe["pipeline_wall_s"].mean())
    per_tile_pipeline_max = float(pipe["pipeline_wall_s"].max())
    sub_per_tile = float(pipe["n_subwindows"].mean())
    rows_per_tile = (float(pipe["rows_state_inventory"].mean())
                     + float(pipe["rows_diplotype"].mean())
                     + float(pipe["rows_phase_pairs"].mean())
                     + float(pipe["rows_repr_info"].mean()))
    workers = 4
    est_wall_h_4w = args.structure_tiles * per_tile_pipeline / workers / 3600.0
    est_wall_h_1w = args.structure_tiles * per_tile_pipeline / 3600.0
    peak_gib = _maxrss_gib()
    est_ram_gib = peak_gib * workers
    est_disk_gib = args.structure_tiles * rows_per_tile * 600.0 / 1073741824.0

    over = est_wall_h_4w > WALL_BUDGET_H or est_ram_gib > 12.0 or est_disk_gib > DISK_BUDGET_GIB
    if over:
        affordable = int(max(1, np.floor(WALL_BUDGET_H * 0.6 * 3600.0 * workers / per_tile_pipeline)))
    else:
        affordable = args.structure_tiles

    report = {
        "run_id": bs.RUN_ID, "protocol_version": bs.PROTOCOL_VERSION,
        "evidence_class": "[실측-신규]", "scope": "[표본] 12 benchmark tiles",
        "stage": "io_benchmark_and_HC-D1_extrapolation",
        "software": {
            "bcftools": subprocess.check_output([vcfio.BCFTOOLS, "--version"]).decode().split("\n")[0],
            "python": sys.version.split()[0], "numpy": np.__version__, "pandas": pd.__version__,
        },
        "benchmark_design": {
            "tiles": int(len(small)), "people": args.bench_people,
            "sample_id_transport": "pipe (/dev/fd), never written to disk",
            "read_command": "bcftools query -r <tile> -S <pipe> -f '%POS\\t%REF\\t%ALT[\\t%GT]\\n'",
            "view_form_also_timed": "bcftools view -H -r <tile> -S <pipe> (reference only)",
            "subwindow_cost_included": True, "markers": 16, "step": 8,
            "outcome_blind": True, "seed": 20260922,
        },
        "per_tile_2048_people": {
            "n_tiles": int(len(small)),
            "median_read_wall_s": float(small["read_wall_s"].median()),
            "median_subwindow_build_wall_s": float(small["subwindow_build_wall_s"].median()),
            "median_total_wall_s": float(small["total_wall_s"].median()),
            "median_primary_variants": float(small["n_primary_variants"].median()),
            "median_records_in_tile": float(small["n_records_in_tile"].median()),
            "median_subwindows": float(small["n_subwindows"].median()),
            "median_calls_per_s": float(small["calls_per_s"].median()),
            "all_separators_pipe": bool(small["sep_all_pipe"].all()),
        },
        "per_tile_full_cohort": {
            "n_people": n_all, "n_tiles_probed": int(len(full)),
            "mean_read_wall_s": float(full["read_wall_s"].mean()),
            "mean_subwindow_build_wall_s": float(full["subwindow_build_wall_s"].mean()),
            "mean_calls_per_s": float(full["calls_per_s"].mean()),
            "read_scaling_vs_2048": float(full["read_wall_s"].mean()
                                          / small["read_wall_s"].head(len(full)).mean()),
            "people_ratio": n_all / float(args.bench_people),
        },
        "per_tile_full_pipeline": {
            "n_tiles_probed": int(len(pipe)),
            "mean_wall_s": per_tile_pipeline, "max_wall_s": per_tile_pipeline_max,
            "mean_subwindows": sub_per_tile,
            "components_mean_s": {k: float(pipe[k].mean()) for k in
                                  ("read_wall_s", "state_wall_s", "contrast_wall_s",
                                   "pair_wall_s", "spectrum_wall_s") if k in pipe.columns},
            "mean_output_rows": rows_per_tile,
        },
        "extrapolation_300_tiles_full_cohort": {
            "workers": workers, "threads_per_worker": 1,
            "est_wall_hours_4_workers": round(est_wall_h_4w, 2),
            "est_wall_hours_1_worker": round(est_wall_h_1w, 2),
            "est_peak_rss_gib_per_worker": round(peak_gib, 2),
            "est_peak_rss_gib_total": round(est_ram_gib, 2),
            "est_new_disk_gib": round(est_disk_gib, 3),
            "genotype_bytes_streamed_gib": round(
                args.structure_tiles * float(full["n_variants_read"].mean()) * n_all * 4.0
                / 1073741824.0, 1),
            "budget_wall_hours": WALL_BUDGET_H, "budget_disk_gib": DISK_BUDGET_GIB,
            "budget_ram_gib": 12.0,
            "within_budget": (not over),
            "status": "OK" if not over else "RESOURCE_LIMITED",
            "proposed_scope_if_limited": int(affordable),
            "scope_rule": ("if over budget: reduce the NUMBER OF TILES only and say so; "
                           "never silently lower people, markers or thresholds"),
        },
        "benchmark_wall_s": round(time.time() - t_start, 1),
        "limitations": [
            "throughput only; no phenotype, no model fit, no performance claim",
            "12 tiles are a [표본]; dense tiles (n_primary>4096) are excluded by the "
            "DENSE_WINDOW_RESOURCE_LIMIT rule and are not represented here",
            "extrapolation assumes 4 concurrent workers and no contention from the "
            "concurrent fset branch; measured load is recorded separately",
        ],
    }
    with open(args.out, "w") as fh:
        json.dump({"report": report,
                   "per_tile_small": rows_small, "per_tile_full": rows_full,
                   "per_tile_view_form": rows_view, "per_tile_pipeline": rows_pipeline},
                  fh, indent=2, ensure_ascii=False, default=float)
    print(json.dumps(report, indent=2, ensure_ascii=False, default=float))
    return 0

if __name__ == "__main__":
    sys.exit(main())
